"""
This contains a range of functions required to correct for the 
Magnetisation Transfer effect visible in the HCP data.
"""

import json
from pathlib import Path
from fsl.wrappers import fslmaths, LOAD, bet, fast
from fsl.data.image import Image
import numpy as np
from .utils import create_dirs, load_json, update_json, setup_logger
from .tissue_masks import generate_tissue_mask
from .distortion_correction import register_fmap, generate_asl_mask
import subprocess
import regtricks as rt
import nibabel as nb

import os
import os.path as op

def generate_asl2struct(asl_vol0, struct, fsdir, reg_dir, verbose):
    """
    Generate the linear transformation between ASL-space and T1w-space
    using FS bbregister. Note that struct is required only for saving 
    the output in the right convention, it is not actually used by 
    bbregister. 
    
    Args:
        asl_vol0: path to first volume of ASL 
        struct: path to T1w image (eg T1w_acdc_restore.nii.gz)
        fsdir: path to subject's FreeSurfer output directory 
        reg_dir: path to registration directory, for output 

    Returns: 
        n/a, file 'asl2struct.mat' will be saved in reg_dir
    """
    # set up logger
    logger_name = "HCPASL.asl2struct_registration"
    log_name = reg_dir/"asl2struct_reg.log"
    logger = setup_logger(logger_name, log_name, "INFO", verbose)

    logger.info("Running generate_asl2struct()")
    logger.info(f"Movable volume: {asl_vol0}")
    logger.info(f"T1w structural image: {struct}")
    logger.info(f"FreeSurfer output directory: {fsdir}")
    logger.info(f"Output directory: {reg_dir}")
    logger.info(f"Verbose: {verbose}")

    # We need to do some hacky stuff to get bbregister to work...
    # Split the path to the FS directory into a fake $SUBJECTS_DIR
    # and subject_id. We temporarily set the environment variable
    # before making the call, and then revert back afterwards  
    new_sd, sid = op.split(fsdir)
    old_sd = os.environ.get('SUBJECTS_DIR')
    pwd = os.getcwd()
    orig_mgz = op.join(fsdir, 'mri', 'orig.mgz')

    # Run inside the regdir. Save the output in fsl format, by default 
    # this targets the orig.mgz, NOT THE T1 IMAGE ITSELF! 
    logger.info(f"Changing directory to {reg_dir}")
    os.chdir(reg_dir)
    omat_path = op.join(reg_dir, "asl2struct.mat")
    logger.info(f"Setting $SUBJECTS_DIR to {new_sd} for bbregister")
    cmd = os.environ['SUBJECTS_DIR'] = new_sd
    cmd = f"$FREESURFER_HOME/bin/bbregister --s {sid} --mov {asl_vol0} --t2 "
    cmd += f"--reg asl2orig_mgz_initial_bbr.dat --fslmat {omat_path} --init-fsl"
    logger.info(f"Running bbregister: {cmd}")
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    while 1:
        retcode = process.poll()
        line = process.stdout.readline().decode("utf-8")
        logger.info(line)
        if line == "" and retcode is not None:
            break
    if retcode != 0:
        logger.info(f"retcode={retcode}")
        logger.exception("Process failed.")

    try:
        asl2orig_fsl = rt.Registration.from_flirt(str(omat_path), str(asl_vol0), str(orig_mgz))
    except RuntimeError as e:
        # final row != [0 0 0 1], round to 5 d.p. and try again
        logger.exception("FSL .mat file has an invalid format. Rounding to 5 d.p.")
        arr = np.loadtxt(omat_path)
        np.savetxt(omat_path, arr, fmt='%.5f')
        asl2orig_fsl = rt.Registration.from_flirt(str(omat_path), str(asl_vol0), str(orig_mgz))

    # Return to original working directory, and flip the FSL matrix to target
    # asl -> T1, not orig.mgz. Save output. 
    os.chdir(pwd)
    if old_sd:
        logger.info(f"Changing $SUBJECTS_DIR back to {old_sd}")
        os.environ['SUBJECTS_DIR'] = old_sd
    logger.info("Converting .mat to target T1w.nii.gz rather than orig.mgz.")
    asl2struct_fsl = asl2orig_fsl.to_flirt(str(asl_vol0), str(struct))
    np.savetxt(op.join(reg_dir, 'asl2struct.mat'), asl2struct_fsl)

def correct_M0(subject_dir, calib_dir, mt_factors, 
               t1w_dir, aslt1w_dir, gradunwarp_dir, topup_dir, 
               wmparc, ribbon, corticallut, subcorticallut, 
               interpolation=3, nobandingcorr=False, 
               outdir="hcp_asl", verbose=False):
    """
    Correct the M0 images.
    
    For each of the subject's two calibration images:
    #. Apply gradient and epi distortion corrections;
    #. Apply MT banding correction;
    #. Estimate registration to structural using FreeSurfer's bbregister;
    #. Use SE-based on the gdc_dc calibration image to obtain the bias-field;
    #. Apply bias correction and MT banding correction to gdc_dc calibration image.
    
    Parameters
    ----------
    subject_dir : pathlib.Path
        Path to the subject's base directory.
    calib_dir : pathlib.Path
        Path to the subject's ASL/Calib directory.
    mt_factors : pathlib.Path
        Path to the empirically estimated MT correction 
        scaling factors.
    t1w_dir : pathlib.Path
        Path to the subject's T1w directory (within the 
        Structural_preproc directory).
    aslt1w_dir : pathlib.Path
        Path to the subject's structural output directory, for  
        example ${SubjectDir}/${OutDir}/ASLT1w.
    gradunwarp_dir : pathlib.Path
        Path to the subject's gradient_unwarp run, for example 
        ${SubjectDir}/${OutDir}/ASL/gradient_unwarp.
    topup_dir : pathlib.Path
        Path to the subject's topup run, for example 
        ${SubjectDir}/${OutDir}/ASL/topup.
    wmparc : pathlib.Path
        Path to the subject's wmparc.nii.gz FreeSurfer output for 
        use in SE-based bias correction.
    ribbon : pathlib.Path
        Path to the subject's ribbon.nii.gz FreeSurfer output for 
        use in SE-based bias correction.
    corticallut : pathlib.Path
        FreeSurferCorticalLabelTableLut.txt for use in SE-based 
        bias correction.
    subcorticallut : pathlib.Path
        FreeSurferSubcorticalLabelTableLut.txt for use in SE-based 
        bias correction.
    interpolation : int, {0, 5}
        Order of interpolation to use when applying transformations.
        Default is 3.
    nobandingcorr : bool, optional
        If this is True, the banding correction options in the 
        pipeline will be switched off. Default is False (i.e. 
        banding corrections are applied by default).
    outdir : str
        Name of the main results directory. Default is 'hcp_asl'.
    """
    logger_name = "HCPASL.correct_M0"
    log_name = subject_dir/outdir/"ASL/Calib/correct_M0.log"
    logger = setup_logger(logger_name, log_name, "INFO", verbose)

    # get calibration image names
    calib0, calib1 = [
        (calib_dir/f"Calib{n}/calib{n}.nii.gz").resolve(strict=True) 
        for n in ("0", "1")
    ]
    
    # get structural image names
    struct_name, struct_brain_name = [
        (t1w_dir/f"T1w_acpc_dc_restore{suf}.nii.gz").resolve(strict=True)
        for suf in ("", "_brain")
    ]

    # generate white matter mask in T1w space for use in registration
    logger.info("Generating white matter mask in T1w space.")
    t1reg_dir = aslt1w_dir/"reg"
    t1reg_dir.mkdir(exist_ok=True, parents=True)
    aparc_aseg = (t1w_dir/"aparc+aseg.nii.gz").resolve(strict=True)
    wmmask_img = generate_tissue_mask(aparc_aseg, "wm")
    wmmask_name = t1reg_dir/"wmmask.nii.gz"
    nb.save(wmmask_img, wmmask_name)

    # load gradient distortion correction warp, fieldmaps and PA epidc warp
    logger.info("Loading gradient and EPI distortion correction warps.")
    gdc_name = (gradunwarp_dir/'fullWarp_abs.nii.gz').resolve(strict=True)
    gdc_warp = rt.NonLinearRegistration.from_fnirt(
        str(gdc_name), str(calib0), str(calib0),
        intensity_correct=True, constrain_jac=(0.01, 100)
    )
    fmap, fmapmag, fmapmagbrain = [topup_dir/f"fmap{ext}.nii.gz" 
                                   for ext in ('', 'mag', 'magbrain')]
    epi_dc_warp = rt.NonLinearRegistration.from_fnirt(coefficients=str(topup_dir/"WarpField_01.nii.gz"),
                                                      src=str(fmap),
                                                      ref=str(fmap),
                                                      intensity_correct=True,
                                                      constrain_jac=(0.01, 100))

    # register fieldmapmag to structural image for use in SE-based later
    logger.info("Getting registration from fmapmag image to structural image.")
    fmap_struct_dir = topup_dir/"fmap_struct_reg"
    Path(fmap_struct_dir).mkdir(exist_ok=True)
    fsdir = (t1w_dir/f"{subject_dir.parts[-1]}_V1_MR").resolve(strict=True)
    generate_asl2struct(fmapmag, struct_name, fsdir, fmap_struct_dir, verbose)
    logger.info("Loading registration from fieldmap to struct.")
    bbr_fmap2struct = rt.Registration.from_flirt(str(fmap_struct_dir/"asl2struct.mat"), 
                                                 src=str(fmapmag), 
                                                 ref=str(struct_name)) 

    # iterate over the two calibration images, applying the corrections to both
    logger.info("Iterating over subject's calibration images.")
    for calib_name in (calib0, calib1):
        # get calib_dir and other info
        calib_dir = calib_name.parent
        calib_name_stem = calib_name.stem.split('.')[0]
        logger.info(f"Processing {calib_name_stem}.")

        # apply gdc and epidc to the calibration image
        logger.info(f"Applying gradient and EPI distortion correction to {calib_name_stem}.")
        gdc_dc_warp = rt.chain(gdc_warp, epi_dc_warp)
        gdc_dc_calib_img = gdc_dc_warp.apply_to_image(str(calib_name), str(calib_name), order=interpolation)
        distcorr_dir = calib_dir/"DistCorr"
        distcorr_dir.mkdir(exist_ok=True)
        gdc_dc_calib_name = distcorr_dir/f"gdc_dc_{calib_name_stem}.nii.gz"
        nb.save(gdc_dc_calib_img, gdc_dc_calib_name)

        # apply mt scaling factors to the gradient distortion-corrected calibration image
        if not nobandingcorr:
            logger.info(f"Applying MT scaling factors to {calib_name_stem}.")
            mt_sfs = np.loadtxt(mt_factors)
            assert (len(mt_sfs) == gdc_dc_calib_img.shape[2])
            mt_gdc_dc_calib_img = nb.nifti1.Nifti1Image(gdc_dc_calib_img.get_fdata()*mt_sfs, 
                                                        gdc_dc_calib_img.affine)
            mtcorr_dir = calib_dir/"MTCorr"
            mtcorr_dir.mkdir(exist_ok=True)
            mt_gdc_dc_calib_name = mtcorr_dir/f"mtcorr_gdc_dc_{calib_name_stem}.nii.gz"
            nb.save(mt_gdc_dc_calib_img, mt_gdc_dc_calib_name)
            calib_corr_name = mt_gdc_dc_calib_name
        else:
            calib_corr_name = gdc_dc_calib_name

        # get registration to structural
        logger.info("Generate registration to structural.")
        generate_asl2struct(calib_corr_name, struct_name, fsdir, distcorr_dir, verbose)
        asl2struct_name = distcorr_dir/"asl2struct.mat"
        asl2struct_reg = rt.Registration.from_flirt(src2ref=str(distcorr_dir/"asl2struct.mat"),
                                                    src=str(calib_corr_name),
                                                    ref=str(struct_name))
        # invert for struct2calib registration
        struct2calib_reg = asl2struct_reg.inverse()
        struct2calib_name = distcorr_dir/"struct2asl.mat"
        np.savetxt(struct2calib_name, struct2calib_reg.to_flirt(str(struct_name), str(calib_corr_name)))
        
        # register fmapmag to calibration image space
        logger.info(f"Registering {fmapmag.stem} to {calib_name_stem}")
        fmap2calib_reg = rt.chain(bbr_fmap2struct, struct2calib_reg)
        fmapmag_calibspc = fmap2calib_reg.apply_to_image(str(fmapmag),
                                                         str(calib_name),
                                                         order=interpolation)
        biascorr_dir = calib_dir/"BiasCorr"
        sebased_dir = biascorr_dir/"SEbased"
        sebased_dir.mkdir(parents=True, exist_ok=True)
        fmapmag_cspc_name = sebased_dir/f"fmapmag_{calib_name_stem}spc.nii.gz"
        nb.save(fmapmag_calibspc, fmapmag_cspc_name)

        # get brain mask in calibration image space
        logger.info("Getting brain mask in calibration image space.")
        fs_brainmask = (t1w_dir/"brainmask_fs.nii.gz").resolve(strict=True)
        aslfs_mask_name = calib_dir/"aslfs_mask.nii.gz"
        aslfs_mask = struct2calib_reg.apply_to_image(src=str(fs_brainmask), 
                                                     ref=str(calib_name),
                                                     order=0)
        aslfs_mask = nb.nifti1.Nifti1Image(np.where(aslfs_mask.get_fdata()>0., 1., 0.),
                                           affine=gdc_dc_calib_img.affine)
        nb.save(aslfs_mask, aslfs_mask_name)

        # get sebased bias estimate
        sebased_cmd = [
            "get_sebased_bias",
            "-i", gdc_dc_calib_name, "-f", fmapmag_cspc_name,
            "-m", aslfs_mask_name, "-o", sebased_dir, 
            "--ribbon", ribbon, "--wmparc", wmparc,
            "--corticallut", corticallut, "--subcorticallut", subcorticallut,
            "--struct2calib", struct2calib_name, "--structural", struct_name, 
            "--debug"
        ]
        logger.info(f"Running SE-based bias estimation on {calib_name_stem}.")
        logger.info(sebased_cmd)
        process = subprocess.Popen(sebased_cmd, stdout=subprocess.PIPE)
        while 1:
            retcode = process.poll()
            line = process.stdout.readline().decode("utf-8")
            logger.info(line)
            if line == "" and retcode is not None:
                break
        if retcode != 0:
            logger.info(f"retcode={retcode}")
            logger.exception("Process failed.")

        # apply dilall to bias estimate
        bias_name = sebased_dir/"sebased_bias_dil.nii.gz"
        dilall_name = biascorr_dir/f"{calib_name_stem}_bias.nii.gz"
        dilall_cmd = ["fslmaths", bias_name, "-dilall", dilall_name]
        subprocess.run(dilall_cmd, check=True)

        # bias correct and mt correct the gdc_dc_calib image
        logger.info(f"Performing bias correction.")
        bias_img = nb.load(dilall_name)
        bc_calib = nb.nifti1.Nifti1Image(gdc_dc_calib_img.get_fdata() / bias_img.get_fdata(),
                                         gdc_dc_calib_img.affine)
        biascorr_name = biascorr_dir / f'{calib_name_stem}_restore.nii.gz'
        nb.save(bc_calib, biascorr_name)

        if not nobandingcorr:
            logger.info(f"Performing MT correction.")
            mt_bc_calib = nb.nifti1.Nifti1Image(bc_calib.get_fdata()*mt_sfs,
                                                bc_calib.affine)
            mtcorr_name = mtcorr_dir / f'{calib_name_stem}_mtcorr.nii.gz'
            nb.save(mt_bc_calib, mtcorr_name)
            calib_corr_name = mtcorr_name
        else:
            calib_corr_name = biascorr_name
        