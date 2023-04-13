"""
Functions for the correction of HCP ASL data.

Corrections to be applied include:
    - Bias-field correction (using the bias field estimated on the 
        subject's calibration image - currently the first 
        calibration image but maybe this should be changed to the 
        mean of the first and the second? and maybe there should be 
        some registration from the calibration image to the mean of 
        the ASL series to register the bias field to the series?)
    - MT correction (correction of the Magnetisation Transfer effect 
        visible in HCP ASL data using empirically estimated scaling 
        coefficients)
    - Saturation recovery (post-MT correction, the saturation 
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
import subprocess
import sys

import nibabel as nb
import numpy as np
import regtricks as rt
from fabber import Fabber, percent_progress
from fsl.wrappers import fslmaths
from fsl.wrappers.flirt import applyxfm, mcflirt
from scipy.ndimage import binary_dilation

from .m0_correction import generate_asl2struct
from .utils import create_dirs, setup_logger

# asl sequence parameters
NTIS = 5
IAF = "tc"
IBF = "tis"
TIS = [1.7, 2.2, 2.7, 3.2, 3.7]
RPTS = [6, 6, 6, 10, 15]
SLICEDT = 0.059
SLICEBAND = 10
NSLICES = 60


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


def _satrecov_worker(control_name, satrecov_dir, tis, rpts, ibf, spatial):
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
    logger_name = "HCPASL.hcp_asl_moco"
    logger = logging.getLogger(logger_name)
    # set options for Fabber run, generic to spatial and non-spatial runs
    options = {
        "data": str(control_name),
        "overwrite": True,
        "noise": "white",
        "ibf": ibf,
        "model": "satrecov",
        "save-mean": True,  # is this needed in nonspatial run?
        "casl": True,
        "slicedt": 0.059,
        "ti1": tis[0],
        "ti2": tis[1],
        "ti3": tis[2],
        "ti4": tis[3],
        "ti5": tis[4],
        "sliceband": 10,
        "bolus": 1.5,
        "rpt1": rpts[0],
        "rpt2": rpts[1],
        "rpt3": rpts[2],
        "rpt4": rpts[3],
        "rpt5": rpts[4],
        "fixa": True,
    }
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
    logger.info(f"Running fabber's satrecov model with options:")
    for key, val in options.items():
        logger.info(f"{key}: {str(val)}")
    # run Fabber
    fab = Fabber()
    run = fab.run(
        options, progress_cb=percent_progress(sys.stdout)
    )  # Basic interaction with the run output
    # info about fabber run
    logger.info("\nOutput data summary")
    for name, data in run.data.items():
        logger.info("%s: %s" % (name, data.shape))
    logger.info("Run finished at: %s" % run.timestamp_str)
    # Write full contents out to a directory
    # load control image to get header for saving
    control_img = nb.load(control_name)
    run.write_to_dir(out_dir, ref_nii=control_img)


def _split_tag_control(asl_name, ntis, iaf, ibf, rpts):
    """
    Split an ASL sequence into its tag and control images.

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
    logger = logging.getLogger("HCPASL.hcp_asl_moco")
    logger.info("Splitting data into tag and control using asl_file.")
    asl_base = asl_name.parent / asl_name.stem.split(".")[0]
    # messy - is there a better way of filling in arguments?
    cmd = [
        "asl_file",
        f"--data={str(asl_name)}",
        f"--ntis={ntis}",
        f"--iaf={iaf}",
        f"--ibf={ibf}",
        "--spairs",
        f"--rpts={rpts[0]},{rpts[1]},{rpts[2]},{rpts[3]},{rpts[4]}",
        f"--out={asl_base}",
    ]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    while 1:
        retcode = process.poll()
        line = process.stdout.readline().decode("utf-8")
        logger.info(line)
        if line == "" and retcode is not None:
            break
    if retcode != 0:
        logger.info(f"retcode={retcode}")
        logger.exception("Process failed.")
    even_name = asl_name.parent / f"{asl_base}_even.nii.gz"
    odd_name = asl_name.parent / f"{asl_base}_odd.nii.gz"
    return even_name, odd_name


def _saturation_recovery(asl_name, results_dir, ntis, iaf, ibf, tis, rpts):
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
    control_name, tag_name = _split_tag_control(asl_name, ntis, iaf, ibf, rpts)
    # satrecov nospatial
    _satrecov_worker(control_name, results_dir, tis, rpts, ibf, spatial=False)
    # satrecov spatial
    _satrecov_worker(control_name, results_dir, tis, rpts, ibf, spatial=True)
    t1_name = results_dir / "spatial/mean_T1t.nii.gz"
    return t1_name


def _fslmaths_med_filter_wrapper(image_name):
    """
    Simple wrapper for fslmaths' median filter function. Applies
    the median filter to `image_name`. Derives and returns the
    name of the filtered image as {image_name}_filt.nii.gz.
    """
    filtered_name = image_name.parent / f'{image_name.stem.split(".")[0]}_filt.nii.gz'
    cmd = ["fslmaths", image_name, "-fmedian", filtered_name]
    subprocess.run(cmd)
    return filtered_name


def _slicetiming_correction(asl_name, t1_name, tis, rpts, slicedt, sliceband, n_slices):
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
    numexp = np.where(t1_data > 0, -tis_array / t1_data, 0)
    num = 1 - np.exp(numexp)
    # divide asl sequence by satrecov model evaluated at actual slice-time
    denexp = np.where(t1_data > 0, -slice_times / t1_data, 0)
    den = 1 - np.exp(denexp)
    # evaluate scaling factors
    stcorr_factors = np.where(den > 0, num / den, 1)
    stcorr_factors_img = nb.nifti1.Nifti1Image(stcorr_factors, affine=asl_img.affine)
    # correct asl series
    stcorr_data = asl_img.get_fdata() * stcorr_factors
    stcorr_img = nb.nifti1.Nifti1Image(stcorr_data, affine=asl_img.affine)
    return stcorr_img, stcorr_factors_img


def _register_param(param_name, transform_dir, reffile, param_reg_name):
    """
    Apply motion estimates to an image to obtain a time series.

    Given a parameter map, `param_name`, and a series of motion
    estimates, `transform_dir`, apply the motion estimates to
    the parameter map and obtain a time series of the map
    in the frame described by the motion estimates.

    Uses fsl's applyxfm.

    Parameters
    ----------
    param_name : pathlib.Path
        Path to the parameter estimate to which we want to apply
        motion.
    transform_dir : pathlib.Path
        Path to a directory containing mcflirt motion estimates.
    reffile : pathlib.Path
        Path for reference image in mcflirt motion estimate.
    param_reg_name : pathlib.Path
        Savename of output file
    """
    # list of transformations in transform_dir
    transforms = sorted(transform_dir.glob("**/*"))
    out_names = []
    for n, transform in enumerate(transforms):
        # naming intermediate file
        if n < 10:
            out_n = (
                param_reg_name.parent
                / f'{param_reg_name.stem.split(".")[0]}_000{n}.nii.gz'
            )
        else:
            out_n = (
                param_reg_name.parent
                / f'{param_reg_name.stem.split(".")[0]}_00{n}.nii.gz'
            )
        out_names.append(out_n)
        applyxfm(str(param_name), str(reffile), str(transform), str(out_n))
        cmd = [
            "applyxfm4D",
            str(param_name),
            str(reffile),
            str(out_n),
            str(transform),
            "-singlematrix",
        ]
        subprocess.run(cmd)
    # merge registered parameter volumes into one time series
    cmd = ["fslmerge", "-t", str(param_reg_name), *out_names]
    subprocess.run(cmd)
    # remove intermediate file names
    for out_n in out_names:
        out_n.unlink()


def single_step_resample_to_asl0(
    subject_dir,
    tis_dir,
    mt_factors,
    bias_name,
    calib_name,
    calib2struct,
    gradunwarp_dir,
    topup_dir,
    t1w_dir,
    cores=1,
    interpolation=3,
    nobandingcorr=False,
    outdir="hcp_asl",
    gd_corr=True,
):
    """
    Full ASL correction and motion estimation pipeline.

    This function performs motion estimation for the HCP ASL data.
    The steps of the pipeline are:
    #. Application of Gradient and EPI distortion correction,
        previously estimated via gradient_unwarp and topup
        respectively;
    #. Bias-field correction using the bias-field estimated from
        the first calibration image;
    #. MT correction using pre-calculated correction factors;
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
    tis_dir : pathlib.Path
        Path to the subject's ASL sequence directory.
    mt_factors : pathlib.Path
        Path to the pre-calculated MT correction scaling factors.
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
        ${SubjectDir}/${OutDir}/ASL/gradient_unwarp.
    topup_dir : pathlib.Path
        Path to the subject's topup run, for example
        ${SubjectDir}/${OutDir}/ASL/topup.
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
    outdir : str
        Name of the main results directory. Default is 'hcp_asl'.
    gd_corr: bool
        Whether to perform gradient distortion correction or not.
        Default is True
    """
    # set up logger
    logger_name = "HCPASL.single_step_resample_to_asl0"
    log_out = subject_dir / outdir / f"ASL/single_step_resample_to_asl0.log"
    logger = setup_logger(logger_name, log_out, "INFO")

    logger.info("Running single_step_resample_to_asl0()")
    logger.info(f"Subject directory: {subject_dir}")
    logger.info(f"ASL TIs directory: {tis_dir}")
    logger.info(f"MT scaling factors: {mt_factors}")
    logger.info(f"SE-based bias name: {bias_name}")
    logger.info(f"Calibration image name: {calib_name}")
    logger.info(f"calib2struct registration: {calib2struct}")
    logger.info(f"gradient_unwarp results directory: {gradunwarp_dir}")
    logger.info(f"Topup results directory: {topup_dir}")
    logger.info(f"T1w Structural Preprocessing directory: {t1w_dir}")
    logger.info(f"Number of CPU cores to use: {cores}")
    logger.info(f"Interpolation order: {interpolation}")
    logger.info(f"Perform banding corrections: {not nobandingcorr}")
    logger.info(f"outdir: {outdir}")

    assert (
        isinstance(cores, int) and cores > 0 and cores <= mp.cpu_count()
    ), f"Number of cores should be an integer from 1-{mp.cpu_count()}."
    assert (
        isinstance(interpolation, int) and interpolation >= 0 and interpolation <= 5
    ), "Order of interpolation should be an integer from 0-5."

    # original ASL series and bias field names
    asl_name = (tis_dir / "tis.nii.gz").resolve(strict=True)

    # create directories for results
    logger.info("Creating results directories.")
    bcorr_dir = tis_dir / "BiasCorr"
    distcorr_dir = tis_dir / "DistCorr/FirstPass"
    moco_dir = tis_dir / "MoCo"
    asln2m0_name = moco_dir / "asln2m0.mat"
    satrecov_dir = tis_dir / "SatRecov"
    create_dirs(
        [tis_dir, bcorr_dir, distcorr_dir, moco_dir, asln2m0_name, satrecov_dir]
    )
    if not nobandingcorr:
        mtcorr_dir = tis_dir / "MTCorr"
        stcorr_dir = tis_dir / "STCorr"
        create_dirs([mtcorr_dir, stcorr_dir])

    # apply gradient distortion correction to the ASL series

    # load gdc warp
    gdc_name = (gradunwarp_dir / "fullWarp_abs.nii.gz").resolve()
    if gd_corr:
        logger.info("Applying gradient distortion correction to original ASL series.")
        gdc_warp = rt.NonLinearRegistration.from_fnirt(
            coefficients=str(gdc_name),
            src=str(asl_name),
            ref=str(asl_name),
            intensity_correct=True,
        )
    # load topup's epi dc warp and fmap2struct registration
    dc_name = (topup_dir / "WarpField_01.nii.gz").resolve(strict=True)
    fmapmag_name = (topup_dir / "fmapmag.nii.gz").resolve(strict=True)
    dc_warp = rt.NonLinearRegistration.from_fnirt(
        coefficients=str(dc_name),
        src=str(fmapmag_name),
        ref=str(fmapmag_name),
        intensity_correct=True,
    )
    fmap2struct_reg = (topup_dir / "fmap_struct_reg/asl2struct.mat").resolve(
        strict=True
    )
    t1w_name = (t1w_dir / "T1w_acpc_dc_restore.nii.gz").resolve(strict=True)
    fmap2struct_reg = rt.Registration.from_flirt(
        src2ref=str(fmap2struct_reg), src=str(fmapmag_name), ref=str(t1w_name)
    )
    # apply gdc to ASL series
    if gd_corr:
        asl_corr = gdc_warp.apply_to_image(
            src=str(asl_name), ref=str(asl_name), order=interpolation, cores=cores
        )
        asl_corr = nb.nifti1.Nifti1Image(
            asl_corr.get_fdata().astype(np.float32), affine=asl_corr.affine
        )
        asl_corr_name = distcorr_dir / "tis_gdc.nii.gz"
        nb.save(asl_corr, asl_corr_name)
    else:
        asl_corr_name = asl_name
        asl_corr = nb.load(asl_name)

    # bias correct the ASL series
    logger.info("Bias-correcting the ASL series.")
    bcorr_img = bcorr_dir / "tis_biascorr.nii.gz"
    fslmaths(str(asl_corr_name)).div(str(bias_name)).run(str(bcorr_img))

    # apply MT scaling factors to the bias-corrected ASL series
    if not nobandingcorr:
        logger.info("MT-correcting the ASL series.")
        mtcorr_name = mtcorr_dir / "tis_mtcorr.nii.gz"
        mt_sfs = np.loadtxt(mt_factors)
        mt_arr = np.repeat(np.tile(mt_sfs, (86, 86, 1))[..., np.newaxis], 86, axis=-1)
        mt_img = nb.nifti1.Nifti1Image(mt_arr, affine=asl_corr.affine)
        biascorr_img = nb.load(bcorr_img)
        assert len(mt_sfs) == biascorr_img.shape[2]
        mtcorr_img = nb.nifti1.Nifti1Image(
            (biascorr_img.get_fdata() * mt_img.get_fdata()).astype(np.float32),
            affine=biascorr_img.affine,
        )
        nb.save(mtcorr_img, mtcorr_name)
        asl_corr = mtcorr_name
    else:
        asl_corr = bcorr_img

    # estimate satrecov model on gradient distortion-, bias- and MT- corrected ASL series
    logger.info("First satrecov model fit.")
    t1_name = _saturation_recovery(asl_corr, satrecov_dir, NTIS, IAF, IBF, TIS, RPTS)
    t1_filt_name = _fslmaths_med_filter_wrapper(t1_name)

    # perform slice-time correction using estimated tissue params
    if not nobandingcorr:
        logger.info("Performing initial ST correction.")
        stcorr_img, stfactors_img = _slicetiming_correction(
            mtcorr_name, t1_filt_name, TIS, RPTS, SLICEDT, SLICEBAND, NSLICES
        )
        stcorr_img, stfactors_img = [
            nb.nifti1.Nifti1Image(img.get_fdata().astype(np.float32), affine=img.affine)
            for img in (stcorr_img, stfactors_img)
        ]
        stcorr_name = stcorr_dir / "tis_stcorr.nii.gz"
        nb.save(stcorr_img, stcorr_name)
        stfactors_name = stcorr_dir / "st_scaling_factors.nii.gz"
        nb.save(stfactors_img, stfactors_name)
        asl_corr = stcorr_name

    # register ASL series to calibration image
    logger.info("Running mcflirt on calibration image and ASL series.")
    reg_name = moco_dir / "initial_registration_TIs.nii.gz"
    mcflirt(
        str(asl_corr), reffile=str(calib_name), mats=True, plots=True, out=str(reg_name)
    )
    # rename mcflirt matrices directory
    orig_mcflirt = moco_dir / "initial_registration_TIs.nii.gz.mat"
    if asln2m0_name.exists():
        shutil.rmtree(asln2m0_name)
    orig_mcflirt.rename(asln2m0_name)

    # get brain mask in ASL0 space
    logger.info("Create brain mask in ASL0 space")
    asln2m0_moco = rt.MotionCorrection.from_mcflirt(
        mats=str(asln2m0_name), src=str(asl_name), ref=str(calib_name)
    )
    fs_brainmask = (t1w_dir / "brainmask_fs.nii.gz").resolve(strict=True)
    calib2struct_reg = rt.Registration.from_flirt(
        src2ref=str(calib2struct), src=str(calib_name), ref=str(fs_brainmask)
    )
    struct2asl0_reg = rt.chain(
        calib2struct_reg.inverse(), asln2m0_moco.transforms[0].inverse()
    )
    aslfs_mask = struct2asl0_reg.apply_to_image(
        src=str(fs_brainmask), ref=str(asl_name), order=0
    )
    aslfs_mask = nb.nifti1.Nifti1Image(
        np.where(aslfs_mask.get_fdata() > 0.0, 1.0, 0.0), affine=aslfs_mask.affine
    )
    aslfs_mask_name = tis_dir / "aslfs_mask.nii.gz"
    nb.save(aslfs_mask, aslfs_mask_name)
    # dilate mask so it isn't so strict and make 4d for application to ASL time series
    mask4d = binary_dilation(
        aslfs_mask.get_fdata()[..., np.newaxis], iterations=2
    ).astype(np.float32)

    # register pre-ST-correction ASLn to ASL0
    logger.info("Apply distortion and motion correction to original ASL series.")
    asl02m0 = asln2m0_moco.transforms[0]
    asln2asl0 = rt.chain(asln2m0_moco, asl02m0.inverse())
    asl02fmap = rt.chain(asl02m0, calib2struct_reg, fmap2struct_reg.inverse())
    dc_asln2asl0 = rt.chain(asln2asl0, asl02fmap, dc_warp, asl02fmap.inverse())
    if gd_corr:
        dc_asln2asl0 = rt.chain(gdc_warp, dc_asln2asl0)
    reg_dc = dc_asln2asl0.apply_to_image(
        str(asl_name), str(calib_name), cores=cores, order=interpolation
    )
    reg_dc = nb.nifti1.Nifti1Image(
        (reg_dc.get_fdata() * mask4d).astype(np.float32), affine=reg_dc.affine
    )
    reg_dc_name = distcorr_dir / "temp_reg_dc_tis.nii.gz"
    nb.save(reg_dc, reg_dc_name)

    # apply moco to MT scaling factors image
    if not nobandingcorr:
        logger.info("Apply motion estimates to the MT scaling factors image")
        mt_reg_img = asln2asl0.apply_to_image(
            src=mt_img, ref=mt_img, cores=cores, order=interpolation
        )
        mt_reg_img = nb.nifti1.Nifti1Image(
            np.where(mask4d != 0.0, mt_reg_img.get_fdata(), 1.0).astype(np.float32),
            affine=mt_reg_img.affine,
        )
        mt_reg_img_name = moco_dir / "reg_mt_scaling_factors.nii.gz"
        nb.save(mt_reg_img, mt_reg_img_name)

    # apply bias-correction to motion- and distortion-corrected ASL series
    logger.info(
        "Apply bias correction to the distortion- and motion-corrected ASL series."
    )
    bias_img = nb.load(bias_name)
    reg_dc_biascorr = nb.nifti1.Nifti1Image(
        ((reg_dc.get_fdata() * mask4d) / bias_img.get_fdata()[..., np.newaxis]).astype(
            np.float32
        ),
        affine=reg_dc.affine,
    )
    reg_dc_biascorr_name = moco_dir / "reg_dc_tis_biascorr.nii.gz"
    nb.save(reg_dc_biascorr, reg_dc_biascorr_name)

    # apply MT-correction
    if not nobandingcorr:
        temp_reg_dc_mtcorr = moco_dir / "temp_reg_dc_tis_mtcorr.nii.gz"
        reg_dc_mtcorr = nb.nifti1.Nifti1Image(
            (reg_dc_biascorr.get_fdata() * mt_reg_img.get_fdata() * mask4d).astype(
                np.float32
            ),
            affine=reg_dc.affine,
        )
        nb.save(reg_dc_mtcorr, temp_reg_dc_mtcorr)
        asl_corr = temp_reg_dc_mtcorr
    else:
        asl_corr = reg_dc_biascorr_name

    # re-estimate satrecov model on distortion- and motion-corrected data
    logger.info("Re-fitting the satrecov model since data has been motion-corrected.")
    satrecov_dir = tis_dir / "SatRecov2"
    stcorr_dir = tis_dir / "STCorr2"
    create_dirs([satrecov_dir, stcorr_dir])
    t1_name = _saturation_recovery(asl_corr, satrecov_dir, NTIS, IAF, IBF, TIS, RPTS)
    t1_filt_name = _fslmaths_med_filter_wrapper(t1_name)

    # register T1t estimates back to the original space of the ASL volumes
    # using current motion estimates for improved motion estimation
    logger.info(
        "Registering the final T1t estimates back to the original space of each ASL volume."
    )
    t1_filt_asln = asln2asl0.inverse().apply_to_image(
        src=str(t1_filt_name), ref=str(asl_name), order=interpolation, cores=cores
    )
    t1_filt_asln = nb.nifti1.Nifti1Image(
        t1_filt_asln.get_fdata().astype(np.float32), affine=t1_filt_asln.affine
    )
    t1_filt_asln_name = t1_filt_name.parent / "mean_T1t_filt_asln.nii.gz"
    nb.save(t1_filt_asln, t1_filt_asln_name)

    # apply gdc and epi dc to the ASL series so it is fully distortion
    # corrected in its original space
    logger.info(
        "Applying distortion corrections to ASL volumes in their original spaces."
    )
    dc_asln2asln = rt.chain(dc_asln2asl0, asln2asl0.inverse())
    dc_asl = dc_asln2asln.apply_to_image(
        src=str(asl_name), ref=str(asl_name), order=interpolation, cores=cores
    )
    dc_asl = nb.nifti1.Nifti1Image(
        dc_asl.get_fdata().astype(np.float32), affine=dc_asl.affine
    )
    dc_asl_name = distcorr_dir / "dc_tis.nii.gz"
    nb.save(dc_asl, dc_asl_name)

    # apply bias correction to the distortion corrected ASL series
    logger.info("Applying bias correction to the distortion corrected ASL series.")
    biascorr_dc_asl = nb.nifti1.Nifti1Image(
        (dc_asl.get_fdata() / bias_img.get_fdata()[..., np.newaxis]).astype(np.float32),
        affine=dc_asl.affine,
    )
    biascorr_dc_asl_name = bcorr_dir / "tis_dc_restore.nii.gz"
    nb.save(biascorr_dc_asl, biascorr_dc_asl_name)

    # reapply banding corrections to  ASL series
    if not nobandingcorr:
        # apply MT correction
        logger.info("Applying MT correction to the distortion corrected ASL series.")
        mtcorr_biascorr_dc_asl = nb.nifti1.Nifti1Image(
            (biascorr_dc_asl.get_fdata() * mt_img.get_fdata()).astype(np.float32),
            affine=biascorr_dc_asl.affine,
        )
        mtcorr_biascorr_dc_asl_name = mtcorr_dir / "tis_dc_restore_mtcorr.nii.gz"
        nb.save(mtcorr_biascorr_dc_asl, mtcorr_biascorr_dc_asl_name)

        # apply refined slice-time correction to distortion corrected ASL series
        logger.info(
            "Apply refined slicetiming correction to the distortion corrected ASL series."
        )
        stcorr_img, stfactors_img = _slicetiming_correction(
            mtcorr_biascorr_dc_asl_name,
            t1_filt_asln_name,
            TIS,
            RPTS,
            SLICEDT,
            SLICEBAND,
            NSLICES,
        )
        stcorr_img, stfactors_img = [
            nb.nifti1.Nifti1Image(img.get_fdata().astype(np.float32), affine=img.affine)
            for img in (stcorr_img, stfactors_img)
        ]
        stcorr_name = stcorr_dir / "tis_dc_restore_mtcorr_stcorr.nii.gz"
        stfactors_name = stcorr_dir / "st_scaling_factors.nii.gz"
        nb.save(stcorr_img, stcorr_name)
        nb.save(stfactors_img, stfactors_name)

        # combined MT and ST scaling factors
        logger.info(
            "Combining the ST and MT scaling factors into one set of scaling factors."
        )
        combined_factors_name = stcorr_dir / "combined_scaling_factors_asln.nii.gz"
        combined_factors_img = nb.nifti1.Nifti1Image(
            (stfactors_img.get_fdata() * mt_img.get_fdata()).astype(np.float32),
            affine=stfactors_img.affine,
        )
        nb.save(combined_factors_img, combined_factors_name)
        asl_corr = stcorr_name
    else:
        combined_factors_img = nb.nifti1.Nifti1Image(
            np.ones_like(biascorr_dc_asl.get_fdata(), dtype=np.float32),
            affine=biascorr_dc_asl.affine,
        )
        combined_factors_name = moco_dir / "combined_scaling_factors_asln.nii.gz"
        nb.save(combined_factors_img, combined_factors_name)
        asl_corr = biascorr_dc_asl_name

    # re-estimate motion correction for distortion and banding corrected ASL series
    logger.info(
        "Re-running mcflirt on calibration image and fully corrected ASL series."
    )
    reg_name = moco_dir / "final_registration_TIs.nii.gz"
    asln2m0_final_name = moco_dir / "asln2m0_final.mat"
    mcflirt(
        str(asl_corr), reffile=str(calib_name), mats=True, plots=True, out=str(reg_name)
    )
    # rename mcflirt matrices directory
    final_mcflirt = moco_dir / "final_registration_TIs.nii.gz.mat"
    if asln2m0_final_name.exists():
        shutil.rmtree(asln2m0_final_name)
    final_mcflirt.rename(asln2m0_final_name)

    # apply gdc, refined moco and epidc to the original ASL series
    logger.info(
        "Applying distortion correction and improved motion estimates to the ASL series."
    )
    asln2m0_final = rt.MotionCorrection.from_mcflirt(
        mats=str(asln2m0_final_name), src=str(asl_name), ref=str(calib_name)
    )
    asl02m0_final = asln2m0_final.transforms[0]
    asln2asl0_final = rt.chain(asln2m0_final, asl02m0_final.inverse())
    asl02fmap_final = rt.chain(
        asl02m0_final, calib2struct_reg, fmap2struct_reg.inverse()
    )
    dc_asln2asl0_final = rt.chain(
        asln2asl0_final, asl02fmap_final, dc_warp, asl02fmap_final.inverse()
    )
    if gd_corr:
        dc_asln2asl0_final = rt.chain(gdc_warp, dc_asln2asl0_final)
    dc_moco_asl_final = dc_asln2asl0_final.apply_to_image(
        src=str(asl_name), ref=str(asl_name), order=interpolation, cores=cores
    )
    dc_moco_asl_final = nb.nifti1.Nifti1Image(
        dc_moco_asl_final.get_fdata().astype(np.float32),
        affine=dc_moco_asl_final.affine,
    )
    gdc_dc_moco_asl_final_name = moco_dir / "tis_dc_moco.nii.gz"
    nb.save(dc_moco_asl_final, gdc_dc_moco_asl_final_name)

    # apply bias correction to the fully distortion and motion corrected ASL series
    logger.info(
        "Applying bias correction to the distortion and motion corrected ASL series."
    )
    dc_moco_asl_final_restore = nb.nifti1.Nifti1Image(
        (dc_moco_asl_final.get_fdata() / bias_img.get_fdata()[..., np.newaxis]).astype(
            np.float32
        ),
        affine=dc_moco_asl_final.affine,
    )
    dc_moco_asl_final_restore_name = tis_dir / "tis_dc_moco_restore.nii.gz"
    nb.save(dc_moco_asl_final_restore, dc_moco_asl_final_restore_name)

    # apply final motion estimates to scaling factors
    logger.info("Applying motion estimates to the scaling factors.")
    combined_factors_moco = asln2asl0_final.apply_to_image(
        src=str(combined_factors_name),
        ref=str(asl_name),
        order=interpolation,
        cores=cores,
    )
    combined_factors_moco = nb.nifti1.Nifti1Image(
        combined_factors_moco.get_fdata().astype(np.float32),
        affine=combined_factors_moco.affine,
    )
    combined_factors_moco_name = tis_dir / "combined_scaling_factors.nii.gz"
    nb.save(combined_factors_moco, combined_factors_moco_name)

    # apply banding corrections to the motion corrected ASL series
    if not nobandingcorr:
        logger.info("Banding correcting the registered ASL series.")
        dc_moco_bandcorr_asl = nb.nifti1.Nifti1Image(
            (
                dc_moco_asl_final_restore.get_fdata()
                * combined_factors_moco.get_fdata()
            ).astype(np.float32),
            affine=dc_moco_asl_final_restore.affine,
        )
        dc_moco_bandcorr_asl_name = tis_dir / "tis_dc_moco_restore_bandcorr.nii.gz"
        nb.save(dc_moco_bandcorr_asl, dc_moco_bandcorr_asl_name)


def single_step_resample_to_aslt1w(
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
    corticallut,
    subcorticallut,
    asl_scaling_factors=None,
    mt_factors=None,
    t1_est=None,
    nobandingcorr=False,
    interpolation=3,
    cores=1.0,
    gd_corr=True,
):
    # set up logger
    logger_name = "HCPASL.single_step_resample_to_aslt1w"
    tis_aslt1w_dir = aslt1w_dir / "TIs"
    tis_aslt1w_dir.mkdir(exist_ok=True)
    log_out = tis_aslt1w_dir / "single_step_resample_to_aslt1w.log"
    logger = setup_logger(logger_name, log_out, "INFO")

    logger.info("Running single_step_resample_to_aslt1w()")
    logger.info(f"ASL series: {asl_name}")
    logger.info(f"Calibration image: {calib_name}")
    logger.info(f"Subject directory: {subject_dir}")
    logger.info(f"T1w directory: {t1w_dir}")
    logger.info(f"Mcflirt motion estimate directory: {moco_dir}")
    logger.info(f"Perfusion image: {perfusion_name}")
    logger.info(f"gradient_unwarp output directory: {gradunwarp_dir}")
    logger.info(f"gd_corr: {gd_corr}")
    logger.info(f"Topup output directory: {topup_dir}")
    logger.info(f"ASLT1w directory: {aslt1w_dir}")
    logger.info(f"ribbon.nii.gz: {ribbon}")
    logger.info(f"wmparc.nii.gz: {wmparc}")
    logger.info(f"Cortical LUT: {corticallut}")
    logger.info(f"Subcortical LUT: {subcorticallut}")
    logger.info(f"ASL scaling factors: {asl_scaling_factors}")
    logger.info(f"MT scaling factors: {mt_factors}")
    logger.info(f"Estimated T1t image: {t1_est}")
    logger.info(f"Perform banding corrections: {not nobandingcorr}")
    logger.info(f"Interpolation order: {interpolation}")
    logger.info(f"Number of CPU cores to use: {cores}")

    # get registration from perfusion image to T1w_acpc_dc_restore
    # using FreeSurfer's bbregister
    logger.info("Getting registration from perfusion image to T1w.")
    reg_dir = tis_aslt1w_dir / "reg"
    reg_dir.mkdir(exist_ok=True, parents=True)
    struct_name = (t1w_dir / "T1w_acpc_dc_restore.nii.gz").resolve(strict=True)
    fsdir = (t1w_dir / subject_dir.stem).resolve(strict=True)
    generate_asl2struct(perfusion_name, struct_name, fsdir, reg_dir)
    asl2struct_reg = rt.Registration.from_flirt(
        src2ref=str(reg_dir / "asl2struct.mat"),
        src=str(perfusion_name),
        ref=str(struct_name),
    )

    # get brain mask in ASL-gridded T1w space
    logger.info("Obtaining brain mask in ASLT1w space.")
    struct_brain_mask = t1w_dir / "brainmask_fs.nii.gz"
    asl_spc = rt.ImageSpace(str(asl_name))
    t1w_spc = rt.ImageSpace(str(struct_brain_mask))
    asl_gridded_t1w_spc = t1w_spc.resize_voxels(asl_spc.vox_size / t1w_spc.vox_size)
    aslt1_brain_mask = rt.Registration.identity().apply_to_image(
        src=str(struct_brain_mask), ref=asl_gridded_t1w_spc, order=0
    )
    aslt1_brain_mask = nb.nifti1.Nifti1Image(
        np.where(aslt1_brain_mask.get_fdata() > 0.0, 1.0, 0.0),
        affine=aslt1_brain_mask.affine,
    )
    aslt1_brain_mask_name = reg_dir / "ASL_grid_T1w_brain_mask.nii.gz"
    nb.save(aslt1_brain_mask, aslt1_brain_mask_name)

    # register a field of view image from ASL0 to ASLT1w space
    perfusion_img = nb.load(perfusion_name)
    fov_asl0 = np.ones_like(perfusion_img.get_fdata())
    fov_aslt1w = asl2struct_reg.apply_to_array(
        data=fov_asl0, src=str(perfusion_name), ref=asl_gridded_t1w_spc, order=0
    )
    fov_aslt1w = nb.nifti1.Nifti1Image(
        np.where(fov_aslt1w > 0.1, 1.0, 0.0), affine=aslt1_brain_mask.affine
    )
    fov_aslt1w_name = reg_dir / "ASL_FoV_mask.nii.gz"
    nb.save(fov_aslt1w, fov_aslt1w_name)

    # use logical_and of the brain mask and FoV mask
    fov_brainmask = nb.nifti1.Nifti1Image(
        np.where(
            np.logical_and(
                fov_aslt1w.get_fdata() > 0, aslt1_brain_mask.get_fdata() > 0
            ),
            1.0,
            0.0,
        ),
        affine=fov_aslt1w.affine,
    )
    fov_brainmask_name = reg_dir / "ASL_FoV_brain_mask.nii.gz"
    nb.save(fov_brainmask, fov_brainmask_name)
    mask4d = fov_brainmask.get_fdata()[..., np.newaxis]

    # register fieldmap magnitude image to ASL-gridded T1w space
    logger.info("Registering fmapmag to ASLT1s space.")
    fmapmag = (topup_dir / "fmapmag.nii.gz").resolve(strict=True)
    fmap2struct_name = (topup_dir / "fmap_struct_reg/asl2struct.mat").resolve(
        strict=True
    )
    fmap2struct_reg = rt.Registration.from_flirt(
        src2ref=str(fmap2struct_name), src=str(fmapmag), ref=str(struct_name)
    )
    fmap_aslt1w = fmap2struct_reg.apply_to_image(
        src=str(fmapmag), ref=str(aslt1_brain_mask_name), order=interpolation
    )
    fmap_aslt1w_name = reg_dir / "fmapmag_aslt1w.nii.gz"
    nb.save(fmap_aslt1w, fmap_aslt1w_name)

    # load gradient_unwarp gradient distortion correction warp
    logger.info("Loading distortion correction warps.")
    gdc_name = (gradunwarp_dir / "fullWarp_abs.nii.gz").resolve()
    if gd_corr:
        logger.info(
            "gradient_unwarp.py was run. Gradient distortion correction will be applied."
        )
        gdc_warp = rt.NonLinearRegistration.from_fnirt(
            coefficients=str(gdc_name),
            src=str(asl_name),
            ref=str(asl_name),
            intensity_correct=True,
        )

    # load topup EPI distortion correction warp
    dc_name = (topup_dir / "WarpField_01.nii.gz").resolve(strict=True)
    dc_warp = rt.NonLinearRegistration.from_fnirt(
        coefficients=str(dc_name),
        src=str(fmapmag),
        ref=str(fmapmag),
        intensity_correct=True,
    )

    # load asl motion correction
    logger.info("Loading motion correction.")
    asln2m0_moco = rt.MotionCorrection.from_mcflirt(
        mats=str(moco_dir), src=str(asl_name), ref=str(calib_name)
    )
    m02asl0 = asln2m0_moco.transforms[0].inverse()
    asln2asl0 = rt.chain(asln2m0_moco, m02asl0)

    # form asl0 and calib0 distortion correction to structural warps
    asl02fmap_reg = rt.chain(asl2struct_reg, fmap2struct_reg.inverse())
    asl0_dc2struct_warp = rt.chain(asl02fmap_reg, dc_warp, fmap2struct_reg)
    asl_moco_dc2struct_warp = rt.chain(asln2asl0, asl0_dc2struct_warp)
    if gd_corr:
        asl_moco_dc2struct_warp = rt.chain(gdc_warp, asl_moco_dc2struct_warp)
    asl_moco2struct = rt.chain(asln2asl0, asl2struct_reg)
    m0_dc2struct_warp = rt.chain(m02asl0, asl0_dc2struct_warp)
    if gd_corr:
        m0_dc2struct_warp = rt.chain(gdc_warp, m0_dc2struct_warp)

    # register calibration image to structural and obtain new
    # SE-based bias estimates

    # register original calibration image to ASL-gridded T1w space
    logger.info("Registering calibration image to ASLT1w space.")
    calib_distcorr_dir = aslt1w_dir / "Calib/Calib0/DistCorr"
    calib_distcorr_dir.mkdir(exist_ok=True, parents=True)
    calib_dc_aslt1w = m0_dc2struct_warp.apply_to_image(
        src=str(calib_name), ref=str(aslt1_brain_mask_name), order=interpolation
    )
    calib_dc_aslt1w_name = calib_distcorr_dir / "calib0_dc.nii.gz"
    nb.save(calib_dc_aslt1w, calib_dc_aslt1w_name)

    # create timing image in calibration image space to register to ASLT1w space
    calib_timing_name = calib_name.parent / "calib_timing.nii.gz"
    create_ti_image(str(calib_name), [8], SLICEBAND, SLICEDT, str(calib_timing_name))

    # register calibration timing image to ASLT1w space
    m02struct = rt.chain(m02asl0, asl2struct_reg)
    calib_aslt1w_timing = m02struct.apply_to_image(
        src=str(calib_timing_name), ref=str(aslt1_brain_mask_name), order=0
    )
    calib_aslt1w_timing_name = aslt1w_dir / "Calib/Calib0/calib_aslt1w_timing.nii.gz"
    nb.save(calib_aslt1w_timing, calib_aslt1w_timing_name)

    # get ASLT1w space SE-based bias estimate
    logger.info("Performing SE-based bias estimation in ASLT1w space.")
    sebased_dir = aslt1w_dir / "Calib/Calib0/SEbased"
    sebased_dir.mkdir(exist_ok=True)
    sebased_cmd = [
        "get_sebased_bias_asl",
        "-i",
        calib_dc_aslt1w_name,
        "-f",
        fmap_aslt1w_name,
        "-m",
        fov_brainmask_name,
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
    bias_name = sebased_dir / "sebased_bias_dil.nii.gz"
    dilall_name = sebased_dir / "sebased_bias_dilall.nii.gz"
    dilall_cmd = ["fslmaths", bias_name, "-dilall", dilall_name]
    subprocess.run(dilall_cmd, check=True)

    # get ASL series in ASL-gridded T1w space along with the scaling factors
    # used to perform banding correction.
    # also apply new bias field estimated above.

    # register ASL series to ASL-gridded T1w space
    logger.info("Registering ASL series to ASLT1w space.")
    distcorr_dir = tis_aslt1w_dir / "DistCorr"
    distcorr_dir.mkdir(exist_ok=True, parents=True)
    asl_dc_moco = asl_moco_dc2struct_warp.apply_to_image(
        src=str(asl_name),
        ref=str(aslt1_brain_mask_name),
        cores=cores,
        order=interpolation,
    )
    asl_dc_moco = nb.nifti1.Nifti1Image(
        (asl_dc_moco.get_fdata() * mask4d).astype(np.float32), affine=asl_dc_moco.affine
    )
    asl_dc_moco_name = distcorr_dir / "tis_dc_moco.nii.gz"
    nb.save(asl_dc_moco, asl_dc_moco_name)

    # register ASL scaling factors to ASL-gridded T1w space
    if asl_scaling_factors:
        logger.info("Registering ASL scaling factors to ASLT1w space.")
        aslt1w_sfs = asl_moco2struct.apply_to_image(
            src=str(asl_scaling_factors),
            ref=str(aslt1_brain_mask_name),
            cores=cores,
            order=interpolation,
        )
        aslt1w_sfs = nb.nifti1.Nifti1Image(
            (aslt1w_sfs.get_fdata() * mask4d).astype(np.float32),
            affine=aslt1w_sfs.affine,
        )
    else:
        aslt1w_sfs = nb.nifti1.Nifti1Image(
            np.ones(asl_dc_moco.shape, dtype=np.float32), affine=asl_dc_moco.affine
        )
    aslt1w_sfs_name = tis_aslt1w_dir / "combined_scaling_factors.nii.gz"
    nb.save(aslt1w_sfs, aslt1w_sfs_name)

    # apply bias field and banding corrections to ASL series
    logger.info(
        "Applying SE-based bias correction and banding corrections to ASL series."
    )
    bias = nb.load(bias_name)
    asl_corr = nb.nifti1.Nifti1Image(
        np.where(
            bias.get_fdata()[..., np.newaxis] != 0,
            (asl_dc_moco.get_fdata() * aslt1w_sfs.get_fdata())
            / bias.get_fdata()[..., np.newaxis],
            0.0,
        ).astype(np.float32),
        affine=asl_dc_moco.affine,
    )
    asl_corr_name = tis_aslt1w_dir / "asl_corr.nii.gz"
    nb.save(asl_corr, asl_corr_name)

    # create TI timing image in ASL space and register to ASL-gridded T1w space
    logger.info("Creating TI image in ASLT1w space for use in oxford_asl.")
    ti_aslt1w_name = tis_aslt1w_dir / "timing_img_aslt1w.nii.gz"
    create_ti_image(str(asl_name), TIS, SLICEBAND, SLICEDT, str(ti_aslt1w_name))
    ti_aslt1w = asl2struct_reg.apply_to_image(
        src=str(ti_aslt1w_name), ref=str(aslt1_brain_mask_name), order=0
    )
    nb.save(ti_aslt1w, ti_aslt1w_name)

    # register T1 image, estimated by the satrecov model, to ASL-gridded T1w space
    logger.info("Registering estimated T1t image to ASLT1w space.")
    t1_est_aslt1w = asl2struct_reg.apply_to_image(
        src=str(t1_est), ref=str(aslt1_brain_mask_name), order=interpolation
    )
    t1_est_aslt1w_name = reg_dir / "mean_T1t_filt_aslt1w.nii.gz"
    nb.save(t1_est_aslt1w, t1_est_aslt1w_name)

    # register calibration image's MT scaling factors to ASL-gridded T1w space
    if mt_factors:
        logger.info(
            "Registering calibration image's MT scaling factors to ASLT1w space."
        )
        mt_sfs = np.loadtxt(mt_factors)
        mt_arr = np.tile(mt_sfs, (86, 86, 1))
        calib_mt_sfs_aslt1w = m02struct.apply_to_array(
            data=mt_arr,
            src=str(calib_name),
            ref=str(aslt1_brain_mask_name),
            order=interpolation,
        )
        calib_mt_sfs_aslt1w = nb.nifti1.Nifti1Image(
            calib_mt_sfs_aslt1w, affine=calib_dc_aslt1w.affine
        )
        calib_mt_sfs_aslt1w_name = (
            calib_distcorr_dir / "calib_mt_scaling_factors.nii.gz"
        )
        nb.save(calib_mt_sfs_aslt1w, calib_mt_sfs_aslt1w_name)

        # perform slicetime correction on the calibration image
        num = 1 - np.exp(
            np.where(t1_est_aslt1w.get_fdata() > 0, -8 / t1_est_aslt1w.get_fdata(), 0)
        )
        den = 1 - np.exp(
            np.where(
                np.logical_and(
                    t1_est_aslt1w.get_fdata() > 0, calib_aslt1w_timing.get_fdata() > 0.1
                ),
                -calib_aslt1w_timing.get_fdata() / t1_est_aslt1w.get_fdata(),
                0,
            )
        )
        calib_aslt1w_stcorr_factors = np.where(den > 0, num / den, 1)
        calib_aslt1w_stcorr_factors_img = nb.nifti1.Nifti1Image(
            calib_aslt1w_stcorr_factors, affine=calib_dc_aslt1w.affine
        )
        calib_aslt1w_stcorr_factors_name = (
            aslt1w_dir / "Calib/Calib0/calib_aslt1w_stcorr_factors.nii.gz"
        )
        nb.save(calib_aslt1w_stcorr_factors_img, calib_aslt1w_stcorr_factors_name)

        # correct the registered, gdc_dc, bias-corrected calibration image for MT effect and ST effect
        calib_biascorr = nb.load(sebased_dir / "calib0_secorr.nii.gz")
        calib_mtcorr = nb.nifti1.Nifti1Image(
            calib_biascorr.get_fdata()
            * calib_mt_sfs_aslt1w.get_fdata()
            * calib_aslt1w_stcorr_factors,
            affine=calib_biascorr.affine,
        )
        calib_mtcorr_name = aslt1w_dir / "Calib/Calib0/calib0_corr_aslt1w.nii.gz"
        nb.save(calib_mtcorr, calib_mtcorr_name)
    else:
        calib_biascorr = nb.load(sebased_dir / "calib0_secorr.nii.gz")
        calib_corr_name = aslt1w_dir / "Calib/Calib0/calib0_corr_aslt1w.nii.gz"
        nb.save(calib_biascorr, calib_corr_name)
