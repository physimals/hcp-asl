"""
This contains a range of functions required to correct for the 
Magnetisation Transfer effect visible in the HCP data.
"""

import logging
from pathlib import Path
import os

import nibabel as nb
import numpy as np
import regtricks as rt

from .tissue_masks import generate_tissue_mask
from .utils import sp_run, ImagePath
from .registration import register_asl2struct


def initial_corrections_calibration(
    subject_dir,
    calib_dir,
    eb_factors,
    t1w_dir,
    aslt1w_dir,
    gradunwarp_dir,
    topup_dir,
    wmparc,
    ribbon,
    interpolation=3,
    nobandingcorr=False,
    gd_corr=True,
):
    """
    For each of the subject's two calibration images:
    #. Apply gradient and susceptibility distortion corrections;
    #. Apply empirical banding banding correction;
    #. Estimate registration to structural using FreeSurfer's bbregister;
    #. Use SE-based on the gdc_sdc calibration image to obtain the bias-field;
    #. Apply bias correction and empirical banding banding correction to gdc_sdc calibration image.

    Parameters
    ----------
    subject_dir : pathlib.Path
        Path to the subject's base directory.
    calib_dir : pathlib.Path
        Path to the subject's ASL/Calib directory.
    eb_factors : pathlib.Path
        Path to the empirically estimated empirical banding correction
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

    hcppipedir = Path(os.environ["HCPPIPEDIR"])
    corticallut = hcppipedir / "global/config/FreeSurferCorticalLabelTableLut.txt"
    subcorticallut = hcppipedir / "global/config/FreeSurferSubcorticalLabelTableLut.txt"

    # get calibration image names
    calib0, calib1 = [
        (calib_dir / f"calib{n}/calib{n}.nii.gz").resolve(strict=True)
        for n in ("0", "1")
    ]

    # get structural image names
    struct_name, struct_brain_name = [
        (t1w_dir / f"T1w_acpc_dc_restore{suf}.nii.gz").resolve(strict=True)
        for suf in ("", "_brain")
    ]

    # generate white matter mask in T1w space for use in registration
    logging.info("Generating white matter mask in T1w space")
    t1reg_dir = aslt1w_dir / "registration"
    t1reg_dir.mkdir(exist_ok=True, parents=True)
    aparc_aseg = (t1w_dir / "aparc+aseg.nii.gz").resolve(strict=True)
    wmmask_img = generate_tissue_mask(aparc_aseg, "wm")
    wmmask_name = t1reg_dir / "wmmask.nii.gz"
    nb.save(wmmask_img, wmmask_name)

    # load gradient distortion correction warp, fieldmaps and PA sdc warp
    logging.info("Loading gradient and susceptibility distortion correction warps")
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
            "gradient_unwarp.py was not run, not applying gradient distortion correction"
        )
    fmap, fmapmag = [topup_dir / f"fmap{ext}.nii.gz" for ext in ("", "mag")]
    sdc_warp = rt.NonLinearRegistration.from_fnirt(
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
    register_asl2struct(fmapmag, struct_name, fsdir, fmap_struct_dir)
    logging.info("Loading registration from fieldmap to struct")
    bbr_fmap2struct = rt.Registration.from_flirt(
        fmap_struct_dir / "asl2struct.mat", src=fmapmag, ref=struct_name
    )

    logging.info("Iterating over subject's calibration images applying corrections")
    for calib_path in (calib0, calib1):
        calib_dir = calib_path.parent
        logging.info(f"Processing {calib_path}")

        # Corrected calibration image - updated as corrections are applied
        calib = ImagePath(calib_path)

        # Do initial corrections for GDC and empirical banding and register to structural
        if gd_corr:
            logging.info(f"Applying gradient distortion correction to {calib.stem}")
            calib_gdc = calib.correct_from_image(
                calib_dir / "gdc",
                "gdc",
                gdc_warp.apply_to_image(calib.path, calib.path, order=interpolation),
            )
        else:
            calib_gdc = calib

        if not nobandingcorr:
            logging.info(
                f"Applying empirical banding scaling factors to {calib_gdc.stem}"
            )
            mt_sfs = np.loadtxt(eb_factors)
            assert len(mt_sfs) == calib_gdc.img.shape[2]
            calib_gdc_eb = calib_gdc.correct_from_data(
                calib_dir / "empirical_banding_correction",
                "eb",
                calib_gdc.img.get_fdata() * mt_sfs,
            )
        else:
            calib_gdc_eb = calib_gdc

        logging.info(f"Getting registration of {calib_gdc_eb.stem} to structural")
        reg_dir = calib_dir / "registration"
        reg_dir.mkdir(exist_ok=True, parents=True)
        register_asl2struct(calib_gdc_eb.path, struct_name, fsdir, reg_dir)
        asl2struct_reg = rt.Registration.from_flirt(
            src2ref=reg_dir / "asl2struct.mat",
            src=calib_gdc_eb.path,
            ref=struct_name,
        )
        struct2calib_reg = asl2struct_reg.inverse()
        struct2calib_name = reg_dir / "struct2asl.mat"
        np.savetxt(
            struct2calib_name,
            struct2calib_reg.to_flirt(struct_name, calib_gdc_eb.path),
        )

        # Now that we have registrations from calib2str and fmap2str, use
        # this to apply gdc and sdc to the original calibration image in a single step
        logging.info(f"Applying combined distortion corrections to {calib.stem}")
        fmap2calib_reg = rt.chain(bbr_fmap2struct, struct2calib_reg)
        sdc_calibspc_warp = rt.chain(fmap2calib_reg.inverse(), sdc_warp, fmap2calib_reg)
        if gd_corr:
            sdc_calibspc_warp = rt.chain(gdc_warp, sdc_calibspc_warp)
            suffix = "gdc_sdc"
        else:
            suffix = "sdc"
        calib_gdc_sdc = calib.correct_from_image(
            calib_dir / "distortion_correction",
            suffix,
            sdc_calibspc_warp.apply_to_image(
                src=calib_path, ref=calib_path, order=interpolation
            ),
        )

        logging.info(
            f"Registering {fmapmag.stem} to {calib.stem} for SE-bae bias estimation"
        )
        fmapmag_calibspc = fmap2calib_reg.apply_to_image(
            fmapmag, calib.path, order=interpolation
        )
        biascorr_dir = calib_dir / "bias_correction"
        sebased_dir = biascorr_dir / "sebased"
        sebased_dir.mkdir(parents=True, exist_ok=True)
        fmapmag_cspc_name = sebased_dir / f"fmapmag_{calib.stem}spc.nii.gz"
        nb.save(fmapmag_calibspc, fmapmag_cspc_name)

        logging.info("Getting brain mask in calibration image space")
        fs_brainmask = (t1w_dir / "brainmask_fs.nii.gz").resolve(strict=True)
        aslfs_mask_name = calib_dir / "asl_fs_mask.nii.gz"
        aslfs_mask = struct2calib_reg.apply_to_image(
            src=fs_brainmask, ref=calib.path, order=1
        )
        aslfs_mask = nb.nifti1.Nifti1Image(
            (aslfs_mask.get_fdata() > 0.5).astype(np.float32), affine=calib.img.affine
        )
        nb.save(aslfs_mask, aslfs_mask_name)

        logging.info(f"Running SE-based bias estimation on {calib_gdc_sdc.stem}")
        sebased_cmd = [
            "get_sebased_bias_asl",
            "-i",
            calib_gdc_sdc.path,
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

        logging.info("Applying dilall to bias correction")
        bias_name = sebased_dir / "sebased_bias_dil.nii.gz"
        bias_final = biascorr_dir / f"{calib.stem}_biasfield.nii.gz"
        dilall_cmd = ["fslmaths", bias_name, "-dilall", bias_final]
        sp_run(dilall_cmd)

        logging.info(f"Performing bias correction on {calib_gdc_sdc.stem}")
        bias_img = nb.load(bias_final)
        calib_gdc_sdc_bc = calib_gdc_sdc.correct_from_data(
            biascorr_dir, "bc", calib_gdc_sdc.img.get_fdata() / bias_img.get_fdata()
        )

        if not nobandingcorr:
            logging.info(
                f"Performing empirical banding correction to {calib_gdc_sdc_bc.stem}"
            )
            calib_gdc_sdc_bc_eb = calib_gdc_sdc_bc.correct_from_data(
                biascorr_dir, "eb", calib_gdc_sdc_bc.img.get_fdata() * mt_sfs
            )
        else:
            calib_gdc_sdc_bc_eb = calib_gdc_sdc_bc

        calib_gdc_sdc_bc_eb.img.to_filename(
            calib_dir / f"{calib.stem}_initial_corrected.nii.gz"
        )
