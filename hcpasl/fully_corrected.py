import logging
from pathlib import Path
import os

import regtricks as rt
import nibabel as nb
import numpy as np

from .utils import SLICEBAND, SLICEDT, TIS, sp_run, ASL_SHAPE, ImagePath
from .registration import register_asl2struct
from .asl_correction import create_ti_image


def fully_correct_asl_calibration_aslt1w(
    asl_name,
    calib_name,
    subject_dir,
    t1w_dir,
    moco_dir,
    perfusion_name,
    gradunwarp_dir,
    topup_dir,
    aslt1w_dir,
    ribbon,
    wmparc,
    asl_scaling_factors=None,
    eb_factors=None,
    t1_est=None,
    interpolation=3,
    nobandingcorr=False,
    cores=1.0,
    gd_corr=True,
):
    """
    Apply all corrections to ASL and calibration images, map into ASLT1w space.

    Parameters
    ----------
    asl_name: pathlib.Path
        Path to the ASL series.
    calib_name: pathlib.Path
        Path to the calibration image.
    subject_dir : pathlib.Path
        Path to the subject's base directory.
    t1w_dir : pathlib.Path
        Path to the subject's T1w dir.
    moco_dir : pathlib.Path
        Path to the subject's motion correction directory
    perfusion_name : pathlib.Path
        Path to the perfusion image in ASL native space.
    gradunwarp_dir : pathlib.Path
        Path to the subject's gradient_unwarp run, for example
        ${SubjectDir}/${OutDir}/ASL/gradient_unwarp.
    topup_dir : pathlib.Path
        Path to the subject's topup run, for example
        ${SubjectDir}/${OutDir}/ASL/topup.
    aslt1w_dir : pathlib.Path
        Path to the subject's ASL T1w directory.
    ribbon : pathlib.Path
        Path to the subject's ribbon.nii.gz FreeSurfer output for
        use in SE-based bias correction.
    wmparc : pathlib.Path
        Path to the subject's wmparc.nii.gz FreeSurfer output for
        use in SE-based bias correction.
    asl_scaling_factors : pathlib.Path
        Path to the ASL scaling factors (optional)
    eb_factors : pathlib.Path
        Path to the pre-calculated empirical banding correction factors.
    t1_est : pathlib.Path
        Path to the T1 tissue estimates.
    interpolation : int, optional
        Order of interpolation to be used by regtricks. This is
        passed to scipy's map_coordinates. See that for more
        information. Default is 3.
    nobandingcorr : bool, optional
        If this is True, the banding correction options in the
        pipeline will be switched off. Default is False (i.e.
        banding corrections are applied by default).
    cores: int
        number of CPU cores to use. Default is 1.
    gd_corr: bool
        Whether to perform gradient distortion correction or not.
        Default is True
    """

    tis_aslt1w_dir = aslt1w_dir / "label_control"
    tis_aslt1w_dir.mkdir(exist_ok=True, parents=True)

    logging.info("Running single_step_resample_to_aslt1w()")
    logging.info(f"ASL series: {asl_name}")
    logging.info(f"Calibration image: {calib_name}")
    logging.info(f"Subject directory: {subject_dir}")
    logging.info(f"T1w directory: {t1w_dir}")
    logging.info(f"Motion estimate directory: {moco_dir}")
    logging.info(f"Perfusion image: {perfusion_name}")
    logging.info(f"gradient_unwarp output directory: {gradunwarp_dir}")
    logging.info(f"gd_corr: {gd_corr}")
    logging.info(f"Topup output directory: {topup_dir}")
    logging.info(f"ASLT1w directory: {aslt1w_dir}")
    logging.info(f"ribbon.nii.gz: {ribbon}")
    logging.info(f"wmparc.nii.gz: {wmparc}")
    logging.info(f"ASL scaling factors: {asl_scaling_factors}")
    logging.info(f"Empirical banding scaling factors: {eb_factors}")
    logging.info(f"Estimated T1t image: {t1_est}")
    logging.info(f"Perform banding corrections: {not nobandingcorr}")
    logging.info(f"Interpolation order: {interpolation}")
    logging.info(f"Number of CPU cores to use: {cores}")

    hcppipedir = Path(os.environ["HCPPIPEDIR"])
    corticallut = hcppipedir / "global/config/FreeSurferCorticalLabelTableLut.txt"
    subcorticallut = hcppipedir / "global/config/FreeSurferSubcorticalLabelTableLut.txt"

    # generate registration from perfusion image to T1w_acpc_dc_restore
    logging.info("Getting registration from perfusion image to T1w.")
    reg_dir = aslt1w_dir / "registration"
    reg_dir.mkdir(exist_ok=True, parents=True)
    struct_name = (t1w_dir / "T1w_acpc_dc_restore.nii.gz").resolve(strict=True)
    fsdir = (t1w_dir / subject_dir.stem).resolve(strict=True)
    register_asl2struct(perfusion_name, struct_name, fsdir, reg_dir)
    asl2struct_reg = rt.Registration.from_flirt(
        src2ref=reg_dir / "asl2struct.mat", src=perfusion_name, ref=struct_name
    )

    # get brain mask in ASL-gridded T1w space
    logging.info("Obtaining brain mask in ASLT1w space.")
    struct_brain_mask = t1w_dir / "brainmask_fs.nii.gz"
    asl_spc = rt.ImageSpace(asl_name)
    t1w_spc = rt.ImageSpace(struct_brain_mask)
    aslt1w_spc = t1w_spc.resize_voxels(asl_spc.vox_size / t1w_spc.vox_size)
    aslt1_brain_mask = rt.Registration.identity().apply_to_image(
        src=struct_brain_mask, ref=aslt1w_spc, order=1
    )
    aslt1_brain_mask = nb.nifti1.Nifti1Image(
        (aslt1_brain_mask.get_fdata() > 0).astype(np.float32),
        affine=aslt1_brain_mask.affine,
    )
    aslt1_brain_mask_name = reg_dir / "ASL_grid_T1w_brain_mask.nii.gz"
    nb.save(aslt1_brain_mask, aslt1_brain_mask_name)

    logging.info("Loading motion correction.")
    asln2calibration_moco = rt.MotionCorrection.from_mcflirt(
        mats=moco_dir, src=asl_spc, ref=calib_name
    )
    m02asl0 = asln2calibration_moco.transforms[0].inverse()
    asln2asl0 = rt.chain(asln2calibration_moco, m02asl0)

    # Load motion-FoV mask prepared earlier and combine with brain mask
    fov_valid_asl = moco_dir.parent / "fov_mask.nii.gz"
    fov_valid_aslt1w = asl2struct_reg.apply_to_image(
        fov_valid_asl, aslt1w_spc, order=1, cores=cores
    )
    fov_valid_aslt1w = fov_valid_aslt1w.get_fdata() > 0.9
    fov_valid_aslt1w_path = reg_dir / "fov_mask.nii.gz"
    aslt1w_spc.save_image(fov_valid_aslt1w, fov_valid_aslt1w_path)
    fov_brain_mask = fov_valid_aslt1w & (aslt1_brain_mask.get_fdata() > 0)
    fov_brainmask_name = reg_dir / "brain_fov_mask.nii.gz"
    aslt1w_spc.save_image(fov_brain_mask, fov_brainmask_name)

    # register fieldmap magnitude image to ASL-gridded T1w space
    logging.info("Registering fmapmag to ASLT1w space.")
    fmapmag = (topup_dir / "fmapmag.nii.gz").resolve(strict=True)
    fmap2struct_name = (topup_dir / "fmap_struct_reg/asl2struct.mat").resolve(
        strict=True
    )
    fmap2struct_reg = rt.Registration.from_flirt(
        src2ref=fmap2struct_name, src=fmapmag, ref=struct_name
    )
    fmap_aslt1w = fmap2struct_reg.apply_to_image(
        src=fmapmag, ref=aslt1w_spc, order=interpolation
    )
    fmap_aslt1w_name = reg_dir / "fmapmag_aslt1w.nii.gz"
    nb.save(fmap_aslt1w, fmap_aslt1w_name)

    # load gradient distortion correction
    logging.info("Loading distortion correction warps.")
    gdc_name = (gradunwarp_dir / "fullWarp_abs.nii.gz").resolve()
    if gd_corr:
        logging.info(
            "gradient_unwarp.py was run. Gradient distortion correction will be applied."
        )
        gdc_warp = rt.NonLinearRegistration.from_fnirt(
            coefficients=gdc_name,
            src=asl_spc,
            ref=asl_spc,
            intensity_correct=True,
        )

    # load susceptibility distortion correction
    sdc_name = (topup_dir / "WarpField_01.nii.gz").resolve(strict=True)
    sdc_warp = rt.NonLinearRegistration.from_fnirt(
        coefficients=sdc_name,
        src=fmapmag,
        ref=fmapmag,
        intensity_correct=True,
    )

    # form asl0 distortion correction to structural warps
    asl02fmap_reg = rt.chain(asl2struct_reg, fmap2struct_reg.inverse())
    asl_moco2struct = rt.chain(asln2asl0, asl2struct_reg)
    asl0_dc2struct_warp = rt.chain(asl02fmap_reg, sdc_warp, fmap2struct_reg)
    asl_moco_dc2struct_warp = rt.chain(asln2asl0, asl0_dc2struct_warp)
    if gd_corr:
        asl_moco_dc2struct_warp = rt.chain(gdc_warp, asl_moco_dc2struct_warp)

    # form calibration distortion correction to structural warps
    calibration_dc2struct_warp = rt.chain(m02asl0, asl0_dc2struct_warp)
    if gd_corr:
        calibration_dc2struct_warp = rt.chain(gdc_warp, calibration_dc2struct_warp)

    # register original calibration image to ASL-gridded T1w space
    logging.info("Registering calibration image to ASLT1w space.")
    calib_name = ImagePath(calib_name)
    calib_dir = aslt1w_dir / "calibration/calib0"
    calib_gdc_sdc = calib_name.correct_from_image(
        calib_dir,
        "gdc_sdc",
        calibration_dc2struct_warp.apply_to_image(
            src=calib_name.path, ref=aslt1w_spc, order=interpolation
        ),
    )

    # Derive new bias estimates on the calibration image in ASLT1w space
    logging.info("Performing SE-based bias estimation in ASLT1w space.")
    sebased_dir = calib_dir / "sebased"
    sebased_dir.mkdir(exist_ok=True, parents=True)
    sebased_cmd = [
        "get_sebased_bias_asl",
        "-i",
        calib_gdc_sdc.path,
        "-f",
        fmap_aslt1w_name,
        "-m",
        aslt1_brain_mask_name,
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
        "--debug",
    ]
    sp_run(sebased_cmd)
    bias_name = sebased_dir / "sebased_bias_dil.nii.gz"
    dilall_name = sebased_dir / "sebased_bias_dilall.nii.gz"
    dilall_cmd = ["fslmaths", bias_name, "-dilall", dilall_name]
    sp_run(dilall_cmd)
    calib_gdc_sdc_bc = ImagePath(sebased_dir / "calib0_secorr.nii.gz")

    # register ASL series to ASL-gridded T1w space
    logging.info("Registering ASL series to ASLT1w space.")
    asl_spc = rt.ImageSpace(asl_name)
    aslt1_spc = rt.ImageSpace(aslt1_brain_mask_name)
    asl_corr_dir = aslt1w_dir / "label_control"
    asl_name = ImagePath(asl_name)
    asl_gdc_mc_sdc = asl_name.correct_from_image(
        asl_corr_dir,
        "gdc_mc_sdc",
        asl_moco_dc2struct_warp.apply_to_image(
            src=asl_name.path,
            ref=aslt1_spc,
            cores=cores,
            order=interpolation,
        ),
    )

    # register ASL scaling factors to ASL-gridded T1w space
    if asl_scaling_factors:
        logging.info("Registering ASL scaling factors to ASLT1w space.")
        asl_sfs = asl_moco2struct.apply_to_image(
            src=asl_scaling_factors,
            ref=aslt1_spc,
            cores=cores,
            order=interpolation,
        )
    else:
        asl_sfs = nb.nifti1.Nifti1Image(
            np.ones(asl_gdc_mc_sdc.img.shape, dtype=np.float32),
            affine=asl_gdc_mc_sdc.img.affine,
        )
    asl_sfs.to_filename(asl_corr_dir / "label_control_scaling_factors.nii.gz")

    # apply bias field and banding corrections to ASL series
    logging.info(
        "Applying SE-based bias correction and banding corrections to ASL series."
    )
    bias = nb.load(bias_name)
    asl_gdc_mc_sdc_bc_st_eb = np.zeros(asl_gdc_mc_sdc.img.shape, dtype=np.float32)
    np.divide(
        asl_gdc_mc_sdc.get_fdata() * asl_sfs.get_fdata(),
        bias.get_fdata()[..., None],
        out=asl_gdc_mc_sdc_bc_st_eb,
        where=(bias.get_fdata()[..., None] != 0),
    )
    asl_gdc_mc_sdc_bc_st_eb = asl_gdc_mc_sdc.correct_from_data(
        asl_corr_dir, "bc", asl_gdc_mc_sdc_bc_st_eb
    )
    asl_gdc_mc_sdc_bc_st_eb.img.to_filename(
        asl_corr_dir / "label_control_corrected.nii.gz"
    )

    # create TI timing image in ASL space and register to ASL-gridded T1w space
    logging.info("Creating TI image in ASLT1w space for use in oxford_asl.")
    ti_aslt1w_name = tis_aslt1w_dir / "timing_img.nii.gz"
    create_ti_image(str(asl_name), TIS, SLICEBAND, SLICEDT, str(ti_aslt1w_name))
    ti_aslt1w = asl2struct_reg.apply_to_image(
        src=ti_aslt1w_name, ref=aslt1_spc, order=1
    )
    nb.save(ti_aslt1w, ti_aslt1w_name)

    # register T1 image, estimated by the satrecov model, to ASL-gridded T1w space
    logging.info("Registering estimated T1t image to ASLT1w space.")
    t1_est_aslt1w = asl2struct_reg.apply_to_image(
        src=t1_est, ref=aslt1_spc, order=interpolation
    )
    t1_est_aslt1w_name = reg_dir / "mean_T1t_filt.nii.gz"
    nb.save(t1_est_aslt1w, t1_est_aslt1w_name)

    # create timing image in calibration space and register to ASL-gridded T1w space
    calib_timing = calib_name.path.parent / "calib_timing.nii.gz"
    create_ti_image(str(calib_name.path), [8], SLICEBAND, SLICEDT, str(calib_timing))
    calib_timing = ImagePath(calib_timing)
    m02struct = rt.chain(m02asl0, asl2struct_reg)
    calib_timing = calib_timing.correct_from_image(
        calib_dir,
        "",
        m02struct.apply_to_image(src=calib_timing.path, ref=aslt1w_spc, order=1),
    )

    # register calibration image's empirical banding scaling factors to ASL-gridded T1w space
    if eb_factors:
        logging.info(
            "Registering calibration image's empirical banding scaling factors to ASLT1w space."
        )
        eb_sfs = np.loadtxt(eb_factors)
        eb_arr = np.tile(eb_sfs, (ASL_SHAPE[0], ASL_SHAPE[1], 1))
        eb_img = m02struct.apply_to_array(
            data=eb_arr,
            src=calib_name.path,
            ref=aslt1_spc,
            order=interpolation,
        )
        eb_img = aslt1_spc.make_nifti(eb_img)
        eb_img.to_filename(calib_dir / "eb_scaling_factors.nii.gz")

        # perform slicetime correction on the calibration image
        num = np.zeros_like(t1_est_aslt1w.get_fdata())
        np.divide(
            -8,
            t1_est_aslt1w.get_fdata(),
            out=num,
            where=(t1_est_aslt1w.get_fdata() > 0),
        )

        num = 1 - np.exp(num)
        den = np.zeros_like(t1_est_aslt1w.get_fdata())
        fltr = (t1_est_aslt1w.get_fdata() > 0) & (calib_timing.get_fdata() > 0.1)
        np.divide(
            -calib_timing.get_fdata(),
            t1_est_aslt1w.get_fdata(),
            out=den,
            where=fltr,
        )
        den = 1 - np.exp(den)
        st_img = np.ones_like(den, dtype=np.float32)
        np.divide(num, den, where=(den > 0), out=st_img)
        st_img = nb.nifti1.Nifti1Image(st_img, affine=calib_gdc_sdc.img.affine)
        st_img.to_filename(calib_dir / "st_scaling_factors.nii.gz")

        # correct the registered, gdc_sdc, bias-corrected calibration image for empirical banding effect and slice-time effect
        calib_gdc_sdc_bc_st_eb = calib_gdc_sdc_bc.correct_from_data(
            calib_dir,
            "st_eb",
            calib_gdc_sdc_bc.get_fdata() * eb_img.get_fdata() * st_img.get_fdata(),
        )
        calib_gdc_sdc_bc_st_eb.img.to_filename(calib_dir / "calib0_corrected.nii.gz")

    else:
        calib_gdc_sdc_bc.img.to_filename(calib_dir / "calib0_corrected.nii.gz")

    # Transform raw ASL and calibration to ASLT1w for QC
    logging.info("Transforming raw ASL to ASLT1w space for QC.")
    asl_uncorr = asl_name.correct_from_image(  # noqa
        asl_corr_dir,
        "uncorrected",
        asl2struct_reg.apply_to_image(
            asl_name.path, aslt1_spc, order=interpolation, cores=cores
        ),
    )

    logging.info("Transforming raw calibration image to ASLT1w space for QC.")
    calib_uncorr = calib_name.correct_from_image(  # noqa
        calib_dir,
        "uncorrected",
        m02struct.apply_to_image(calib_name.path, ref=aslt1_spc, order=interpolation),
    )
