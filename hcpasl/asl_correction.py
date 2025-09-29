"""
Functions for the correction of HCP ASL data.

Corrections to be applied include:
    - Bias-field correction (using the bias field estimated on the
        subject's calibration image - currently the first
        calibration image but maybe this should be changed to the
        mean of the first and the second? and maybe there should be
        some registration from the calibration image to the mean of
        the ASL series to register the bias field to the series?)
    - Empirical banding correction (correction of the Magnetisation Transfer effect
        visible in HCP ASL data using empirically estimated scaling
        coefficients)
    - Saturation recovery (post-empirical banding correction, the saturation
        recovery expected in ASL is much more visible, giving a
        residual banding effect)
    - Initial slice-timing correction (done using the parameter
        estimates obtained above. Initial because the satrecov
        model is initially estimated upon data which hasn't been
        motion corrected - we account for this later)
    - Motion estimation
    - Second slice-timing correction (use the motion estimates
        from above to align the parameter estimates with individual
        frames of the ASL series prior to correction. Perform
        slice-timing correction again using these aligned
        parameter maps.)
    - Registration (apply motion estimates to the final slice-timing
        corrected ASL series from the above step to obtain final,
        fully-corrected, registered ASL series)
"""

import logging
import multiprocessing as mp
import shutil
import sys

import nibabel as nb
import numpy as np
import regtricks as rt
from fabber import Fabber, percent_progress
from fsl.wrappers import fslmaths
from fsl.wrappers.flirt import mcflirt
from scipy.ndimage import binary_dilation

from hcpasl.utils import (
    ImagePath,
    sp_run,
    make_motion_fov_mask,
    AslParams,
)


def create_ti_image(asl, tis, sliceband, slicedt, outname, repeats=None):
    """
    Create a 4D series of actual TIs at each voxel.

    Args:
        asl: path to image in the space we wish to create the TI series
        tis: list of TIs in the acquisition
        sliceband: number of slices per band in the acquisition
        slicedt: time taken to acquire each slice
        outname: path to which the ti image is saved
        repeats: list of repeats for each TI. If not provided, assume 1 volume

    Returns:
        n/a, file outname is created in output directory
    """

    asl_spc = rt.ImageSpace(asl)
    n_slice = asl_spc.size[2]
    slice_in_band = np.tile(np.arange(0, sliceband), n_slice // sliceband).reshape(
        1, 1, n_slice, 1
    )
    ti_array = np.array([np.tile(x, asl_spc.size) for x in tis]).transpose(1, 2, 3, 0)
    ti_array = ti_array + (slice_in_band * slicedt)

    if repeats:
        ti_array = np.concatenate(
            [
                np.repeat(ti_array[..., ti, np.newaxis], 2 * r, axis=-1)
                for ti, r in enumerate(repeats)
            ],
            axis=-1,
        )

    rt.ImageSpace.save_like(asl, ti_array, outname)


def _satrecov_worker(control_name, satrecov_dir, params: AslParams, spatial):
    """
    Wrapper for fabber's saturation recovery model.

    Parameters
    ----------
    control_name : pathlib.Path
        Path to the control images.
    satrecov_dir : pathlib.Path
        Parent directory for the satrecov results. Results from
        this will be stored either in {`satrecov_dir`}/spatial
        or {`satrecov_dir`}/nospatial depending on value of
        'spatial`.
    tis : list
        TIs for the ASL sequence.
    rpts : list
        Number of repeats for each TI in the sequence.
    ibf : str
        Input format of the sequence.
    spatial : bool
        Choose whether to run fabber in spatial (True) or
        non-spatial (False) mode.
    """
    # set options for Fabber run, generic to spatial and non-spatial runs
    options = {
        "data": str(control_name),
        "overwrite": True,
        "noise": "white",
        "ibf": params.ibf,
        "model": "satrecov",
        "save-mean": True,  # is this needed in nonspatial run?
        "casl": True,
        "slicedt": params.slicedt,
        "sliceband": params.sliceband,
        "bolus": params.bolus,
        "fixa": True,
    }

    # variable number of TIs and RPTs
    options["ntis"] = len(params.tis)
    for i, ti in enumerate(params.tis, start=1):
        options[f"ti{i}"] = ti
    for i, rp in enumerate(params.rpts, start=1):
        options[f"rpt{i}"] = rp

    # spatial or non-spatial specific options
    spatial_dir = satrecov_dir / "spatial"
    nospatial_dir = satrecov_dir / "nospatial"
    if spatial:
        out_dir = str(spatial_dir)
        extra_options = {
            "method": "spatialvb",
            "output": out_dir,
            "continue-from-mvn": str(nospatial_dir / "finalMVN.nii.gz"),
        }
    else:
        out_dir = str(nospatial_dir)
        extra_options = {"method": "vb", "output": out_dir, "save-mvn": True}
    options.update(extra_options)
    logging.info("Running fabber's satrecov model with options:")
    for key, val in options.items():
        logging.info(f"{key}: {str(val)}")
    # run Fabber
    fab = Fabber()
    run = fab.run(
        options, progress_cb=percent_progress(sys.stdout)
    )  # Basic interaction with the run output
    # info about fabber run
    logging.info("\nOutput data summary")
    for name, data in run.data.items():
        logging.info("%s: %s" % (name, data.shape))
    logging.info("Run finished at: %s" % run.timestamp_str)
    # Write full contents out to a directory
    # load control image to get header for saving
    control_img = nb.load(control_name)
    run.write_to_dir(out_dir, ref_nii=control_img)


def split_asl_label_control(asl_name, ntis, iaf, ibf, rpts):
    """
    Split an ASL sequence into its label and control images.

    Parameters
    ----------
    asl_name : pathlib.Path
        Path to the ASL series to be split.
    ntis : int
        Number of TIs in the ASL sequence.
    iaf : str
        ("tc", "ct", "diff")
    ibf : str
        ("tis", "rpts")
    rpts : list
        List of number of repeats at each TI.

    Returns
    -------
    even_name : pathlib.Path
    odd_name : pathlib.Path
    """
    logging.info("Splitting data into label and control using asl_file.")
    asl_base = asl_name.parent / asl_name.stem.split(".")[0]
    # messy - is there a better way of filling in arguments?
    rpts_str = ",".join([str(r) for r in rpts])
    cmd = [
        "asl_file",
        f"--data={str(asl_name)}",
        f"--ntis={ntis}",
        "--iaf=tc",
        f"--ibf={ibf}",
        "--spairs",
        f"--rpts={rpts_str}",
        f"--out={asl_base}",
    ]
    sp_run(cmd)
    even_name = asl_name.parent / f"{asl_base}_even.nii.gz"
    odd_name = asl_name.parent / f"{asl_base}_odd.nii.gz"
    return even_name, odd_name


def fit_satrecov_model(asl_name, results_dir, params: AslParams):
    """
    Use Fabber's `satrecov` model to estimate a T1 map.

    Split the ASL sequence into tag and control images. Fit the
    `satrecov` model on the control images, first with Fabber's
    spatial mode off, then with it on.

    Parameters
    ----------
    asl_name : pathlib.Path
        Path for the ASL series on which the model will be
        estimated.
    results_dir : pathlib.Path
        Directory in which to save the `nospatial` and `spatial`
        parameter estimates.
    """
    # obtain control images of ASL series
    control_name, tag_name = split_asl_label_control(
        asl_name, params.ntis, "tc", params.ibf, params.rpts
    )
    # satrecov nospatial
    _satrecov_worker(control_name, results_dir, params, spatial=False)
    # satrecov spatial
    _satrecov_worker(control_name, results_dir, params, spatial=True)
    t1_name = results_dir / "spatial/mean_T1t.nii.gz"
    return t1_name


def fslmaths_median_filter(image_name):
    """
    Simple wrapper for fslmaths' median filter function. Applies
    the median filter to `image_name`. Derives and returns the
    name of the filtered image as {image_name}_filt.nii.gz.
    """
    filtered_name = image_name.parent / f'{image_name.stem.split(".")[0]}_filt.nii.gz'
    cmd = ["fslmaths", image_name, "-fmedian", filtered_name]
    sp_run(cmd)
    return filtered_name


def apply_slicetime_correction(
    asl_name, t1_name, tis, rpts, slicedt, sliceband, n_slices
):
    """
    Correct an ASL sequence for the slice-timing effect.

    satrecov_model: S(t) = M_0 * (1 - exp{-t/T1t})

    Divides the value in a voxel by the above model evaluated at
    t = TI + n*slicedt where n is the slice number of the voxel
    within a band. This accounts for the voxel being imaged at a
    different time than its supposed TI.

    Multiplies the value in a voxel by the above model evaluated
    at t = TI, i.e. scales the values as if they had been imaged
    at the TI that was specified in the ASL sequence.

    Parameters
    ----------
    asl_name : pathlib.Path
        Path to the ASL sequence.
    t1_name : pathlib.Path
        Path to the T1 estimate. This T1 estimate can be a single
        volume or a time-series. If it is a time-series, the
        number of time points must match that of the ASL series.
    tis : list
        List of TIs used of the sequence.
    rpts : list
        List of repeats. This is the number of repeats for each
        TI.
    slicedt : float
        Time to acquire a slice in the sequence.
    sliceband : int
        Number of slices in each band in the multi-band
        acquisition.
    n_slices : int
        Number of slices in the sequence.

    Returns
    -------
    stcorr_img : nb.nifti1.Nifti1Image
        Slice-time corrected ASL sequence.
    stcorr_factors_img : nb.nifti1.Nifti1Image
        Scaling factors used to perform the slice-time
        correction.
    """
    # timing information for scan
    # supposed measurement time of slice
    tis_array = np.repeat(np.array(tis), 2 * np.array(rpts)).reshape(1, 1, 1, -1)

    # actual measurment time of slice
    slice_numbers = np.tile(np.arange(0, sliceband), n_slices // sliceband).reshape(
        1, 1, -1, 1
    )
    slice_times = tis_array + (slicedt * slice_numbers)

    # load images
    asl_img = nb.load(asl_name)
    t1_img = nb.load(t1_name)

    # check dimensions of t1 image to see if time series or not
    if t1_img.ndim == 3:
        t1_data = t1_img.get_fdata()[..., np.newaxis]
    elif t1_img.ndim == 4:
        t1_data = t1_img.get_fdata()

    # multiply asl sequence by satrecov model evaluated at TI
    numexp = np.zeros(asl_img.shape)
    np.divide(-tis_array, t1_data, out=numexp, where=(t1_data > 0))
    num = 1 - np.exp(numexp)

    # divide asl sequence by satrecov model evaluated at actual slice-time
    denexp = np.zeros(asl_img.shape)
    np.divide(-slice_times, t1_data, out=denexp, where=(t1_data > 0))
    den = 1 - np.exp(denexp)

    # evaluate scaling factors
    stcorr_factors = np.ones_like(num, dtype=np.float32)
    np.divide(num, den, out=stcorr_factors, where=(den > 0))
    stcorr_factors_img = nb.nifti1.Nifti1Image(stcorr_factors, affine=asl_img.affine)

    # correct asl series
    stcorr_data = asl_img.get_fdata() * stcorr_factors
    stcorr_img = nb.nifti1.Nifti1Image(stcorr_data, affine=asl_img.affine)
    return stcorr_img, stcorr_factors_img


def initial_corrections_asl(
    subject_dir,
    label_control_dir,
    eb_factors,
    bias_name,
    calib_name,
    calib2struct,
    gradunwarp_dir,
    topup_dir,
    t1w_dir,
    cores=1,
    interpolation=3,
    nobandingcorr=False,
    gd_corr=True,
    params: AslParams = None,
):
    """
    Full ASL correction and motion estimation pipeline.

    This function performs motion estimation for the HCP ASL data.
    The steps of the pipeline are:
    #. Application of gradient and susceptibility distortion correction,
        previously estimated via gradient_unwarp and topup
        respectively;
    #. Bias-field correction using the bias-field estimated from
        the first calibration image;
    #. Empirical banding correction using pre-calculated correction factors;
    #. An initial fit of the `satrecov` model on the bias and
    MT-corrected series;
    #. Median filtering of the estimated T1 map;
    #. Initial slice-time correction using the median filtered
    T1 map;
    #. Motion estimation;
    #. Apply the motion estimates to the pre-slicetime-correction
    ASL series;
    #. Second fit of the 'satrecov' model on motion-corrected
    series;
    #. Refined slice-timing correction;
    #. Use calib2struct and asl2m0 registrations to get brain mask
    in ASL0 space for use in oxford_asl;
    #. Create timing image in ASL0 space for use with oxford_asl.

    Parameters
    ----------
    subject_dir : pathlib.Path
        Path to the subject's base directory.
    label_control_dir : pathlib.Path
        Path to the subject's ASL sequence directory.
    eb_factors : pathlib.Path
        Path to the pre-calculated empirical banding correction scaling factors.
    bias_name : str
        Path to the SE-based bias field, estimated on the
        calibration image.
    calib_name : str
        Path to the corrected calibration image to which we will
        register using mcflirt.
    calib2struct : str
        Path to the .mat giving the registration from the
        calibration image to the T1w structural image.
    gradunwarp_dir : pathlib.Path
        Path to the subject's gradient_unwarp run, for example
        ${SubjectDir}/ASL/gradient_unwarp.
    topup_dir : pathlib.Path
        Path to the subject's topup run, for example
        ${SubjectDir}/ASL/topup.
    cores : int, optional
        Number of cores regtricks will use. Default is 1.
    interpolation : int, optional
        Order of interpolation to be used by regtricks. This is
        passed to scipy's map_coordinates. See that for more
        information. Default is 3.
    nobandingcorr : bool, optional
        If this is True, the banding correction options in the
        pipeline will be switched off. Default is False (i.e.
        banding corrections are applied by default).
    gd_corr: bool
        Whether to perform gradient distortion correction or not.
        Default is True
    """
    logging.info("Running single_step_resample_to_asl0()")
    logging.info(f"Subject directory: {subject_dir}")
    logging.info(f"ASL label-control directory: {label_control_dir}")
    logging.info(f"empirical banding scaling factors: {eb_factors}")
    logging.info(f"SE-based bias name: {bias_name}")
    logging.info(f"Calibration image name: {calib_name}")
    logging.info(f"calib2struct registration: {calib2struct}")
    logging.info(f"gradient_unwarp results directory: {gradunwarp_dir}")
    logging.info(f"Topup results directory: {topup_dir}")
    logging.info(f"T1w Structural Preprocessing directory: {t1w_dir}")
    logging.info(f"Number of CPU cores to use: {cores}")
    logging.info(f"Interpolation order: {interpolation}")
    logging.info(f"Perform banding corrections: {not nobandingcorr}")

    assert (
        isinstance(cores, int) and cores > 0 and cores <= mp.cpu_count()
    ), f"Number of cores should be an integer from 1-{mp.cpu_count()}."
    assert (
        isinstance(interpolation, int) and interpolation >= 0 and interpolation <= 5
    ), "Order of interpolation should be an integer from 0-5."

    # original ASL series and bias field names
    asl_name = ImagePath(label_control_dir / "label_control.nii.gz")

    # create directories for results
    logging.info("Creating results directories.")
    bcorr_dir = label_control_dir / "bias_correction"
    sdc_dir = label_control_dir / "susceptibility_distortion_correction/first"
    moco_dir = label_control_dir / "motion_correction"
    asln2calibration_name = moco_dir / "asln2m0.mat"
    satrecov_dir = label_control_dir / "saturation_recovery/first"
    for d in [
        label_control_dir,
        bcorr_dir,
        sdc_dir,
        moco_dir,
        asln2calibration_name,
        satrecov_dir,
    ]:
        d.mkdir(parents=True, exist_ok=True)
    if not nobandingcorr:
        eb_dir = label_control_dir / "empirical_banding_correction"
        stcorr_dir = label_control_dir / "slicetime_correction/first"
        for d in [eb_dir, stcorr_dir]:
            d.mkdir(parents=True, exist_ok=True)

    # apply gradient distortion correction to the ASL series
    # load gdc warp
    gdc_name = (gradunwarp_dir / "fullWarp_abs.nii.gz").resolve()
    asl_spc = rt.ImageSpace(str(asl_name))
    if gd_corr:
        logging.info("Applying gradient distortion correction to original ASL series.")
        gdc_warp = rt.NonLinearRegistration.from_fnirt(
            coefficients=gdc_name,
            src=asl_spc,
            ref=asl_spc,
            intensity_correct=True,
        )
    # load topup's sdc warp and fmap2struct registration
    dc_name = (topup_dir / "WarpField_01.nii.gz").resolve(strict=True)
    fmapmag_name = (topup_dir / "fmapmag.nii.gz").resolve(strict=True)
    sdc_warp = rt.NonLinearRegistration.from_fnirt(
        coefficients=dc_name,
        src=fmapmag_name,
        ref=fmapmag_name,
        intensity_correct=True,
    )
    fmap2struct_reg = (topup_dir / "fmap_struct_reg/asl2struct.mat").resolve(
        strict=True
    )
    t1w_name = (t1w_dir / "T1w_acpc_dc_restore.nii.gz").resolve(strict=True)
    fmap2struct_reg = rt.Registration.from_flirt(
        src2ref=fmap2struct_reg, src=fmapmag_name, ref=t1w_name
    )
    # apply gdc to ASL series
    if gd_corr:
        asl_gdc = gdc_warp.apply_to_image(
            src=asl_name.path, ref=asl_spc, order=interpolation, cores=cores
        )
        asl_gdc = asl_name.correct_from_image(sdc_dir, "gdc", asl_gdc)
    else:
        asl_gdc = asl_name

    # bias correct the ASL series
    logging.info("Bias-correcting the ASL series.")
    asl_gdc_bc = bcorr_dir / f"{asl_gdc.stem}_bc.nii.gz"
    fslmaths(str(asl_gdc)).div(str(bias_name)).run(str(asl_gdc_bc))
    asl_gdc_bc = ImagePath(asl_gdc_bc)

    # apply empirical banding scaling factors to the bias-corrected ASL series
    if not nobandingcorr:
        logging.info("Applying empirical banding correction to the ASL series.")
        eb_sfs = np.loadtxt(eb_factors)
        # tile scaling factors in x,y and along slice (z) dimension = len(eb_sfs)
        eb_arr = np.tile(eb_sfs, (asl_spc.size[0], asl_spc.size[1], 1))
        eb_img = asl_spc.make_nifti(eb_arr)
        eb_img.to_filename(eb_dir / "eb_scaling_factors.nii.gz")
        eb_img = ImagePath(eb_dir / "eb_scaling_factors.nii.gz")
        assert len(eb_sfs) == asl_gdc_bc.img.shape[2]
        asl_gdc_bc_eb = asl_gdc_bc.correct_from_data(
            eb_dir, "eb", asl_gdc_bc.img.get_fdata() * eb_arr[..., None]
        )
    else:
        asl_gdc_bc_eb = asl_gdc_bc

    # estimate satrecov model on gradient distortion-, bias- and MT- corrected ASL series
    logging.info("First satrecov model fit.")
    if params is None:
        raise ValueError("AslParams must be provided to initial_corrections_asl")
    t1_name = fit_satrecov_model(asl_gdc_bc_eb.path, satrecov_dir, params)
    t1_filt_name = fslmaths_median_filter(t1_name)

    # perform slice-time correction using estimated tissue params
    if not nobandingcorr:
        logging.info("Performing initial slice-time correction.")
        stcorr_img, stfactors_img = apply_slicetime_correction(
            asl_gdc_bc_eb.path,
            t1_filt_name,
            params.tis,
            params.rpts,
            params.slicedt,
            params.sliceband,
            params.nslices,
        )
        asl_gdc_bc_eb_st = asl_gdc_bc_eb.correct_from_image(
            stcorr_dir, "st", stcorr_img
        )
        stfactors_img.to_filename(stcorr_dir / "st_scaling_factors.nii.gz")

    else:
        asl_gdc_bc_eb_st = asl_gdc_bc_eb

    # register ASL series to calibration image
    logging.info("Running mcflirt on calibration image and ASL series.")
    mc_img = moco_dir / "first/label_control_mc.nii.gz"
    mcflirt(
        str(asl_gdc_bc_eb_st.path),
        reffile=str(calib_name),
        mats=True,
        plots=True,
        out=str(mc_img),
        stages=4,
    )
    # rename mcflirt matrices directory and load transform
    orig_mcflirt = mc_img.with_suffix(mc_img.suffix + ".mat")
    if asln2calibration_name.exists():
        shutil.rmtree(asln2calibration_name)
    orig_mcflirt.rename(asln2calibration_name)

    # Generate motion-FoV mask in ASL0 space
    asl_spc = rt.ImageSpace(asl_gdc_bc_eb.path)
    calib_spc = rt.ImageSpace(calib_name)
    asln2calibration_moco = rt.MotionCorrection.from_mcflirt(
        mats=asln2calibration_name,
        src=asl_spc,
        ref=calib_spc,
    )
    asl02m0 = asln2calibration_moco.transforms[0]
    asln2asl0 = rt.chain(asln2calibration_moco, asl02m0.inverse())
    fov_mask_asl = make_motion_fov_mask(asln2asl0, asl_spc, asl_spc, cores=cores)
    fov_mask_asl_path = moco_dir / "fov_mask_initial.nii.gz"
    nb.save(fov_mask_asl, fov_mask_asl_path)

    # get brain mask in ASL0 space
    # dilate so not too strict
    logging.info("Create brain mask in ASL0 space")
    fs_brainmask = (t1w_dir / "brainmask_fs.nii.gz").resolve(strict=True)
    calib2struct_reg = rt.Registration.from_flirt(
        src2ref=calib2struct, src=calib_name, ref=fs_brainmask
    )
    struct2asl0_reg = rt.chain(
        calib2struct_reg.inverse(), asln2calibration_moco.transforms[0].inverse()
    )
    aslfs_mask = struct2asl0_reg.apply_to_image(src=fs_brainmask, ref=asl_spc, order=1)
    asl_fs_brainmask = label_control_dir / "brain_mask.nii.gz"
    nb.save(aslfs_mask, asl_fs_brainmask)
    aslfs_mask = binary_dilation(aslfs_mask.get_fdata(), iterations=1).astype(
        np.float32
    )

    # Combine the FS and FoV masks
    # make 4d for application to ASL time series
    asl_mask = (aslfs_mask > 0) & (fov_mask_asl.get_fdata() > 0)
    asl_mask = asl_spc.make_nifti(asl_mask)
    asl_mask_name = label_control_dir / "brain_fov_mask_initial.nii.gz"
    nb.save(asl_mask, asl_mask_name)
    mask4d = asl_mask.get_fdata()[..., None]

    # Start afresh with the raw ASL series, and apply the motion correction
    # and susceptibility distortion correction to get ASLn aligned with ASL0
    logging.info(
        "Apply susceptibility distortion and motion correction to original ASL series."
    )
    asl02fmap = rt.chain(asl02m0, calib2struct_reg, fmap2struct_reg.inverse())
    sdc_asln2asl0 = rt.chain(asln2asl0, asl02fmap, sdc_warp, asl02fmap.inverse())
    if gd_corr:
        dc_asln2asl0 = rt.chain(gdc_warp, sdc_asln2asl0)
    else:
        dc_asln2asl0 = sdc_asln2asl0

    asl_mc_sdc = asl_name.correct_from_image(
        sdc_dir,
        "sdc",
        dc_asln2asl0.apply_to_image(
            asl_name.path, calib_name, cores=cores, order=interpolation
        ),
    )

    # apply moco to empirical banding scaling factors image
    if not nobandingcorr:
        logging.info(
            "Apply motion estimates to the empirical banding scaling factors image"
        )
        eb_mc = asln2asl0.apply_to_array(
            eb_arr, src=eb_img.img, ref=eb_img.img, cores=cores, order=interpolation
        )
        eb_mc = np.where(mask4d != 0.0, eb_mc, 1.0).astype(np.float32)
        eb_mc = eb_img.correct_from_data(moco_dir, "mc", eb_mc)

    # apply bias-correction to motion- and distortion-corrected ASL series
    logging.info(
        "Apply bias correction to the susceptibility distortion and motion corrected ASL series."
    )
    bias_img = nb.load(bias_name)
    asl_mc_sdc_bc = asl_mc_sdc.correct_from_data(
        moco_dir, "bc", asl_mc_sdc.get_fdata() / bias_img.get_fdata()[..., None]
    )

    # apply empirical banding correction
    if not nobandingcorr:
        asl_mc_sdc_bc_eb = asl_mc_sdc_bc.correct_from_data(
            moco_dir, "eb", asl_mc_sdc_bc.get_fdata() * eb_mc.get_fdata()
        )
    else:
        asl_mc_sdc_bc_eb = asl_mc_sdc_bc

    # re-estimate satrecov model on distortion- and motion-corrected data
    logging.info("Re-fitting the satrecov model since data has been motion-corrected.")
    satrecov_dir = label_control_dir / "saturation_recovery/second"
    stcorr_dir = label_control_dir / "slicetime_correction/second"
    for d in [satrecov_dir, stcorr_dir]:
        d.mkdir(parents=True, exist_ok=True)
    t1_name = fit_satrecov_model(asl_mc_sdc_bc_eb.path, satrecov_dir, params)
    t1_filt_name = fslmaths_median_filter(t1_name)

    # register T1t estimates back to the original space of the ASL volumes
    # using current motion estimates for improved motion estimation
    logging.info(
        "Registering the final T1t estimates back to the original space of each ASL volume."
    )
    t1_filt_asln = asln2asl0.inverse().apply_to_image(
        src=t1_filt_name, ref=asl_spc, order=interpolation, cores=cores
    )
    t1_filt_asln = nb.nifti1.Nifti1Image(
        t1_filt_asln.get_fdata().astype(np.float32), affine=t1_filt_asln.affine
    )
    t1_filt_asln_name = t1_filt_name.parent / "mean_T1t_filt_asln.nii.gz"
    nb.save(t1_filt_asln, t1_filt_asln_name)

    # apply sdc to the ASL series so it is fully distortion corrected in its original space
    logging.info(
        "Applying susceptibility distortion correction the original ASL series"
    )
    sdc_asln2asln = rt.chain(sdc_asln2asl0, asln2asl0.inverse())
    asl_sdc = asl_name.correct_from_image(
        sdc_dir,
        "sdc",
        sdc_asln2asln.apply_to_image(
            src=asl_name.path, ref=asl_spc, order=interpolation, cores=cores
        ),
    )

    # apply bias correction to the distortion corrected ASL series
    logging.info(
        "Applying bias correction to the susceptibility distortion corrected ASL series."
    )
    asl_sdc_bc = asl_sdc.correct_from_data(
        bcorr_dir, "bc", (asl_sdc.get_fdata() / bias_img.get_fdata()[..., None])
    )

    # Reapply banding corrections to ASL series
    if not nobandingcorr:
        logging.info(
            "Applying empirical banding correction to the distortion corrected ASL series."
        )
        asl_sdc_bc_eb = asl_sdc_bc.correct_from_data(
            eb_dir, "eb", (asl_sdc_bc.get_fdata() * eb_img.get_fdata()[..., None])
        )

        # apply refined slice-time correction to distortion corrected ASL series
        logging.info(
            "Apply refined slicetiming correction to the distortion corrected ASL series."
        )
        stcorr_img, stfactors_img = apply_slicetime_correction(
            asl_sdc_bc_eb.path,
            t1_filt_asln_name,
            params.tis,
            params.rpts,
            params.slicedt,
            params.sliceband,
            params.nslices,
        )
        asl_sdc_bc_eb_st = asl_sdc_bc_eb.correct_from_image(
            stcorr_dir, "st", stcorr_img
        )
        stfactors_name = stcorr_dir / "st_scaling_factors.nii.gz"
        stfactors_img.to_filename(stfactors_name)

        # combined empirical banding and slice-time scaling factors
        logging.info(
            "Combining the slice-time and empirical banding scaling factors into one set of scaling factors."
        )
        combined_factors_name = stcorr_dir / "combined_scaling_factors_asln.nii.gz"
        combined_factors_img = nb.nifti1.Nifti1Image(
            (stfactors_img.get_fdata() * eb_img.get_fdata()[..., None]).astype(
                np.float32
            ),
            affine=stfactors_img.affine,
        )
        nb.save(combined_factors_img, combined_factors_name)

    else:
        logging.info(
            "No banding correction performed; combined scaling factors will be ones."
        )
        combined_factors_img = nb.nifti1.Nifti1Image(
            np.ones_like(asl_sdc_bc.get_fdata(), dtype=np.float32),
            affine=asl_sdc_bc.img.affine,
        )
        combined_factors_name = moco_dir / "combined_scaling_factors_asln.nii.gz"
        nb.save(combined_factors_img, combined_factors_name)
        asl_sdc_bc_eb_st = asl_sdc_bc

    # re-estimate motion correction for distortion and banding corrected ASL series
    logging.info("Re-running mcflirt on calibration image and corrected ASL series.")
    mc_out = moco_dir / "second/label_control_mc.nii.gz"
    asln2calibration_final_name = moco_dir / "asln2calibration_final.mat"
    mcflirt(
        str(asl_sdc_bc_eb_st.path),
        reffile=str(calib_name),
        mats=True,
        plots=True,
        out=str(mc_out),
        stages=4,
    )
    # rename mcflirt matrices directory
    final_mcflirt = mc_out.with_suffix(mc_out.suffix + ".mat")
    if asln2calibration_final_name.exists():
        shutil.rmtree(asln2calibration_final_name)
    final_mcflirt.rename(asln2calibration_final_name)

    # Update the motion FoV mask
    asl_spc = rt.ImageSpace(asl_sdc_bc_eb_st.path)
    calib_spc = rt.ImageSpace(calib_name)
    asln2calibration_final = rt.MotionCorrection.from_mcflirt(
        asln2calibration_final_name, src=asl_spc, ref=calib_spc
    )
    asl02calibration_final = asln2calibration_moco.transforms[0]
    asln2asl0_final = rt.chain(asln2calibration_moco, asl02m0.inverse())
    fov_mask_asl = make_motion_fov_mask(asln2asl0, asl_spc, asl_spc, cores=cores)
    fov_mask_asl_path = moco_dir / "fov_mask.nii.gz"
    nb.save(fov_mask_asl, fov_mask_asl_path)

    # Combine again with the FS brain mask
    asl_mask = (aslfs_mask > 0) & (fov_mask_asl.get_fdata() > 0)
    asl_mask = asl_spc.make_nifti(asl_mask)
    asl_mask_name = label_control_dir / "brain_fov_mask.nii.gz"
    nb.save(asl_mask, asl_mask_name)
    mask4d = asl_mask.get_fdata()[..., None]

    # apply gdc, refined moco and sdc to the original ASL series
    logging.info(
        "Applying distortion correction and improved motion estimates to the original ASL series."
    )
    asln2calibration_final = rt.MotionCorrection.from_mcflirt(
        mats=asln2calibration_final_name, src=asl_spc, ref=calib_name
    )
    asl02calibration_final = asln2calibration_final.transforms[0]
    asln2asl0_final = rt.chain(asln2calibration_final, asl02calibration_final.inverse())
    asl02fmap_final = rt.chain(
        asl02calibration_final, calib2struct_reg, fmap2struct_reg.inverse()
    )
    dc_asln2asl0_final = rt.chain(
        asln2asl0_final, asl02fmap_final, sdc_warp, asl02fmap_final.inverse()
    )
    if gd_corr:
        dc_asln2asl0_final = rt.chain(gdc_warp, dc_asln2asl0_final)

    asl_gdc_mc_sdc = asl_name.correct_from_image(
        label_control_dir,
        "gdc_mc_sdc",
        dc_asln2asl0_final.apply_to_image(
            src=asl_name.path, ref=asl_spc, order=interpolation, cores=cores
        ),
    )

    # apply bias correction to the fully distortion and motion corrected ASL series
    logging.info(
        "Applying bias correction to the distortion and motion corrected ASL series."
    )
    asl_gdc_mc_sdc_bc = asl_gdc_mc_sdc.correct_from_data(
        label_control_dir,
        "bc",
        asl_gdc_mc_sdc.get_fdata() / bias_img.get_fdata()[..., None],
    )

    # apply final motion estimates to scaling factors
    logging.info("Applying motion estimates to the scaling factors.")
    combined_factors_moco = asln2asl0_final.apply_to_image(
        src=combined_factors_name,
        ref=asl_spc,
        order=interpolation,
        cores=cores,
    )
    combined_factors_moco.to_filename(
        label_control_dir / "combined_scaling_factors_mc.nii.gz"
    )

    # apply banding corrections to the motion corrected ASL series
    if not nobandingcorr:
        logging.info("Banding correcting the registered ASL series.")
        asl_gdc_mc_sdc_bc_st_eb = asl_gdc_mc_sdc_bc.correct_from_data(  # noqa
            label_control_dir,
            "eb_st",
            asl_gdc_mc_sdc_bc.get_fdata() * combined_factors_moco.get_fdata(),
        )
        asl_gdc_mc_sdc_bc_st_eb.img.to_filename(
            label_control_dir / "label_control_corrected.nii.gz"
        )
    else:
        asl_gdc_mc_sdc_bc.img.to_filename(
            label_control_dir / "label_control_corrected.nii.gz"
        )
