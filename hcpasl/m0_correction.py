"""
This contains a range of functions required to correct for the 
Magnetisation Transfer effect visible in the HCP data.
"""
import logging
import os.path as op
from pathlib import Path

import nibabel as nb
import numpy as np
import regtricks as rt

from .tissue_masks import generate_tissue_mask
from .utils import sp_run, ImagePath

def generate_asl2struct(asl_vol0, struct, fsdir, reg_dir):
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

    logging.info("Running generate_asl2struct()")
    logging.info(f"Movable volume: {asl_vol0}")
    logging.info(f"T1w structural image: {struct}")
    logging.info(f"FreeSurfer output directory: {fsdir}")
    logging.info(f"Output directory: {reg_dir}")

    # We need to do some hacky stuff to get bbregister to work...
    # Split the path to the FS directory into a fake $SUBJECTS_DIR
    # and subject_id. We temporarily set the environment variable
    # before making the call, and then revert back afterwards
    new_sd, sid = op.split(fsdir)
    orig_mgz = op.join(fsdir, "mri", "orig.mgz")

    # Save the output in fsl format, by default
    # this targets the orig.mgz, NOT THE T1 IMAGE ITSELF!
    logging.info(f"Running bbregister in {reg_dir}; setting $SUBJECTS_DIR to {new_sd}")
    omat_path = op.join(reg_dir, "asl2struct.mat")
    cmd = f"$FREESURFER_HOME/bin/bbregister --s {sid} --mov {asl_vol0} --t2 "
    cmd += f"--reg asl2orig_mgz_initial_bbr.dat --fslmat {omat_path} --init-fsl"
    fslog_name = op.join(reg_dir, "asl2orig_mgz_initial_bbr.dat.log")
    logging.info(f"FreeSurfer's bbregister log: {fslog_name}")
    sp_run(cmd, shell=True, env={"SUBJECTS_DIR": new_sd}, cwd=reg_dir)

    # log final .dat transform
    with open(omat_path, "r") as f:
        lines = f.readlines()
        for line in lines:
            logging.info(line)

    # log minimum registration cost
    mincost = np.loadtxt(op.join(reg_dir, "asl2orig_mgz_initial_bbr.dat.mincost"))
    logging.info(f"bbregister's mincost: {mincost[0]:4f}")

    # convert .dat to .mat
    try:
        asl2orig_fsl = rt.Registration.from_flirt(omat_path, asl_vol0, orig_mgz)
    except RuntimeError as e:
        # final row != [0 0 0 1], round to 5 d.p. and try again
        logging.warning("FSL .mat file has an invalid format. Rounding to 5 d.p.")
        arr = np.loadtxt(omat_path)
        np.savetxt(omat_path, arr, fmt="%.5f")
        asl2orig_fsl = rt.Registration.from_flirt(omat_path, asl_vol0, orig_mgz)

    # Return to original working directory, and flip the FSL matrix to target
    # asl -> T1, not orig.mgz. Save output.
    logging.info("Converting .mat to target T1w.nii.gz rather than orig.mgz")
    asl2struct_fsl = asl2orig_fsl.to_flirt(asl_vol0, struct)
    np.savetxt(op.join(reg_dir, "asl2struct.mat"), asl2struct_fsl)


def correct_M0(
    subject_dir,
    calib_dir,
    mt_factors,
    t1w_dir,
    aslt1w_dir,
    gradunwarp_dir,
    topup_dir,
    wmparc,
    ribbon,
    corticallut,
    subcorticallut,
    interpolation=3,
    nobandingcorr=False,
    outdir="hcp_asl",
    gd_corr=True,
):
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
        example ${SubjectDir}/${OutDir}/T1w/ASL.
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
    gd_corr: bool
        Whether to perform gradient distortion correction or not.
        Default is True
    """

    # get calibration image names
    calib0, calib1 = [
        (calib_dir / f"Calib{n}/calib{n}.nii.gz").resolve(strict=True)
        for n in ("0", "1")
    ]

    # get structural image names
    struct_name, struct_brain_name = [
        (t1w_dir / f"T1w_acpc_dc_restore{suf}.nii.gz").resolve(strict=True)
        for suf in ("", "_brain")
    ]

    # generate white matter mask in T1w space for use in registration
    logging.info("Generating white matter mask in T1w space")
    t1reg_dir = aslt1w_dir / "reg"
    t1reg_dir.mkdir(exist_ok=True, parents=True)
    aparc_aseg = (t1w_dir / "aparc+aseg.nii.gz").resolve(strict=True)
    wmmask_img = generate_tissue_mask(aparc_aseg, "wm")
    wmmask_name = t1reg_dir / "wmmask.nii.gz"
    nb.save(wmmask_img, wmmask_name)

    # load gradient distortion correction warp, fieldmaps and PA epidc warp
    logging.info("Loading gradient and EPI distortion correction warps")
    if gd_corr:
        gdc_name = (gradunwarp_dir / "fullWarp_abs.nii.gz").resolve()
        logging.info(f"gradient_unwarp.py was run, loading {gdc_name}")
        gdc_warp = rt.NonLinearRegistration.from_fnirt(
            coefficients=gdc_name,
            src=calib0,
            ref=calib0,
            intensity_correct=True,
        )
    else:
        logging.info(
            f"gradient_unwarp.py was not run, not applying gradient distortion correction"
        )
    fmap, fmapmag, fmapmagbrain = [
        topup_dir / f"fmap{ext}.nii.gz" for ext in ("", "mag", "magbrain")
    ]
    epidc_warp = rt.NonLinearRegistration.from_fnirt(
        coefficients=topup_dir / "WarpField_01.nii.gz",
        src=fmap,
        ref=fmap,
        intensity_correct=True,
    )

    # register fieldmapmag to structural image for use in SE-based later
    logging.info("Getting registration from fmapmag image to structural image")
    fmap_struct_dir = topup_dir / "fmap_struct_reg"
    Path(fmap_struct_dir).mkdir(exist_ok=True, parents=True)
    fsdir = (t1w_dir / subject_dir.stem).resolve(strict=True)
    generate_asl2struct(fmapmag, struct_name, fsdir, fmap_struct_dir)
    logging.info("Loading registration from fieldmap to struct")
    bbr_fmap2struct = rt.Registration.from_flirt(
        fmap_struct_dir / "asl2struct.mat", src=fmapmag, ref=struct_name
    )

    logging.info("Iterating over subject's calibration images applying corrections")
    for calib_path in (calib0, calib1):
        calib_dir = calib_path.parent
        calib = ImagePath(calib_path)
        logging.info(f"Processing {calib}")

        # Corrected calibration image - updated as corrections are applied
        calib_corr = calib

        # Do initial corrections for GDC and MT and register to structural
        if gd_corr:
            logging.info(f"Applying gradient distortion correction to {calib_corr.stem}")
            calib_corr = calib_corr.correct_from_image(
                calib_dir / "GDC", "gdc", 
                gdc_warp.apply_to_image(calib_path, calib_path, order=interpolation)
            )

        if not nobandingcorr:
            logging.info(f"Applying MT scaling factors to {calib_corr.stem}")
            mt_sfs = np.loadtxt(mt_factors)
            assert len(mt_sfs) == calib_corr.img.shape[2]
            calib_corr = calib_corr.correct_from_data(
                calib_dir / "MTCorr", "mtcorr",
                calib_corr.img.get_fdata() * mt_sfs
            )

        logging.info(f"Getting registration of {calib_corr.stem} to structural")
        reg_dir = calib_dir / "StrucReg"
        reg_dir.mkdir(exist_ok=True, parents=True)
        generate_asl2struct(calib_corr.path, struct_name, fsdir, reg_dir)
        asl2struct_reg = rt.Registration.from_flirt(
            src2ref=reg_dir / "asl2struct.mat",
            src=calib_corr.path,
            ref=struct_name,
        )
        struct2calib_reg = asl2struct_reg.inverse()
        struct2calib_name = reg_dir / "struct2asl.mat"
        np.savetxt(
            struct2calib_name,
            struct2calib_reg.to_flirt(struct_name, calib_corr.path),
        )

        # Now that we have registrations from calib2str and fmap2str, use
        # this to apply gdc and epidc to the original calibration image in a single step
        logging.info(f"Applying combined distortion corrections to {calib.stem}")
        fmap2calib_reg = rt.chain(bbr_fmap2struct, struct2calib_reg)
        dc_calibspc_warp = rt.chain(
            fmap2calib_reg.inverse(), epidc_warp, fmap2calib_reg
        )
        if gd_corr:
            dc_calibspc_warp = rt.chain(gdc_warp, dc_calibspc_warp)
            suffix = "gdc_edc"
        else:
            suffix = "edc"
        calib_corr = calib.correct_from_image(
            calib_dir / "DistCorr", suffix,
            dc_calibspc_warp.apply_to_image(src=calib_path, ref=calib_path, order=interpolation)
        )

        logging.info(f"Registering {fmapmag.stem} to {calib.stem} for SEbased bias estimation")
        fmapmag_calibspc = fmap2calib_reg.apply_to_image(
            fmapmag, calib.path, order=interpolation
        )
        biascorr_dir = calib_dir / "BiasCorr"
        sebased_dir = biascorr_dir / "SEbased"
        sebased_dir.mkdir(parents=True, exist_ok=True)
        fmapmag_cspc_name = sebased_dir / f"fmapmag_{calib.stem}spc.nii.gz"
        nb.save(fmapmag_calibspc, fmapmag_cspc_name)

        logging.info("Getting brain mask in calibration image space")
        fs_brainmask = (t1w_dir / "brainmask_fs.nii.gz").resolve(strict=True)
        aslfs_mask_name = calib_dir / "aslfs_mask.nii.gz"
        aslfs_mask = struct2calib_reg.apply_to_image(
            src=fs_brainmask, ref=calib.path, order=1
        )
        aslfs_mask = nb.nifti1.Nifti1Image(
            (aslfs_mask.get_fdata() > 0.5).astype(np.float32), affine=calib.img.affine
        )
        nb.save(aslfs_mask, aslfs_mask_name)

        logging.info(f"Running SE-based bias estimation on {calib_corr.stem}")
        sebased_cmd = [
            "get_sebased_bias_asl",
            "-i",
            calib_corr.path,
            "-f",
            fmapmag_cspc_name,
            "-m",
            aslfs_mask_name,
            "-o",
            sebased_dir,
            "--ribbon",
            ribbon,
            "--wmparc",
            wmparc,
            "--corticallut",
            corticallut,
            "--subcorticallut",
            subcorticallut,
            "--struct2calib",
            struct2calib_name,
            "--structural",
            struct_name,
            "--debug",
        ]
        sp_run(sebased_cmd)

        logging.info(f"Applying dilall to bias correction")
        bias_name = sebased_dir / "sebased_bias_dil.nii.gz"
        dilall_name = biascorr_dir / f"{calib.stem}_bias.nii.gz"
        dilall_cmd = ["fslmaths", bias_name, "-dilall", dilall_name]
        sp_run(dilall_cmd)

        logging.info(f"Performing bias correction on {calib_corr.stem}")
        bias_img = nb.load(dilall_name)
        calib_corr = calib_corr.correct_from_data(
            biascorr_dir, "bc", 
            calib_corr.img.get_fdata() / bias_img.get_fdata()
        )

        if not nobandingcorr:
            logging.info(f"Performing MT correction to {calib_corr.stem}")
            calib_corr = calib_corr.correct_from_data(
                biascorr_dir, "mt", 
                calib_corr.img.get_fdata() * mt_sfs
            )
