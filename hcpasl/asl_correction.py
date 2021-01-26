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

from .utils import create_dirs
from .m0_mt_correction import load_json, update_json
from .distortion_correction import generate_asl_mask
from fsl.wrappers import fslmaths, LOAD
from fsl.wrappers.flirt import mcflirt, applyxfm, applyxfm4D
from fsl.data.image import Image
import nibabel as nb
from fabber import Fabber, percent_progress
import sys
from pathlib import Path
import shutil
import subprocess
import numpy as np
import regtricks as rt
import multiprocessing as mp
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
    # set options for Fabber run, generic to spatial and non-spatial runs
    options = {
        'data': str(control_name),
        'overwrite': True,
        'noise': 'white',
        'ibf': ibf,
        'model': 'satrecov',
        'save-mean': True, # is this needed in nonspatial run?
        'casl': True,
        'slicedt': 0.059,
        'ti1': tis[0],
        'ti2': tis[1],
        'ti3': tis[2],
        'ti4': tis[3],
        'ti5': tis[4],
        'sliceband': 10,
        'bolus': 1.5,
        'rpt1': rpts[0],
        'rpt2': rpts[1],
        'rpt3': rpts[2],
        'rpt4': rpts[3],
        'rpt5': rpts[4],
        'fixa': True
    }
    # spatial or non-spatial specific options
    spatial_dir = satrecov_dir / 'spatial'
    nospatial_dir = satrecov_dir / 'nospatial'
    if spatial:
        out_dir = str(spatial_dir)
        extra_options = {
            'method': 'spatialvb',
            'output': out_dir,
            'continue-from-mvn': str(nospatial_dir / 'finalMVN.nii.gz')
        }
    else:
        out_dir = str(nospatial_dir)
        extra_options = {
            'method': 'vb',
            'output': out_dir,
            'save-mvn': True
        }
    options.update(extra_options)
    # run Fabber
    fab = Fabber()
    run = fab.run(options, progress_cb=percent_progress(sys.stdout))# Basic interaction with the run output
    # info about fabber run
    print("\nOutput data summary")
    for name, data in run.data.items():
        print("%s: %s" % (name, data.shape))
    print("Run finished at: %s" % run.timestamp_str)
    # Write full contents out to a directory
    # load control image to get header for saving
    control_img = Image(str(control_name))
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
    asl_base = asl_name.parent / asl_name.stem.split('.')[0]
    # messy - is there a better way of filling in arguments?
    cmd = [
        'asl_file',
        f'--data={str(asl_name)}',
        f'--ntis={ntis}',
        f'--iaf={iaf}',
        f'--ibf={ibf}',
        '--spairs',
        f'--rpts={rpts[0]},{rpts[1]},{rpts[2]},{rpts[3]},{rpts[4]}',
        f'--out={asl_base}'
    ]
    subprocess.run(cmd)
    even_name = asl_name.parent / f'{asl_base}_even.nii.gz'
    odd_name = asl_name.parent / f'{asl_base}_odd.nii.gz'
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
    t1_name = results_dir / 'spatial/mean_T1t.nii.gz'
    return t1_name

def _fslmaths_med_filter_wrapper(image_name):
    """
    Simple wrapper for fslmaths' median filter function. Applies 
    the median filter to `image_name`. Derives and returns the 
    name of the filtered image as {image_name}_filt.nii.gz.
    """
    filtered_name = image_name.parent / f'{image_name.stem.split(".")[0]}_filt.nii.gz'
    cmd = [
        'fslmaths',
        image_name,
        '-fmedian',
        filtered_name
    ]
    subprocess.run(cmd)
    return filtered_name

def _slicetiming_correction(
    asl_name, t1_name, tis, rpts, 
    slicedt, sliceband, n_slices
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
    stcorr_img : fsl.image.Image
        Slice-time corrected ASL sequence.
    stcorr_factors_img : fsl.image.Image
        Scaling factors used to perform the slice-time 
        correction.
    """
    # timing information for scan
    # supposed measurement time of slice
    tis_array = np.repeat(np.array(tis), 2*np.array(rpts)).reshape(1, 1, 1, -1)
    # actual measurment time of slice
    slice_numbers = np.tile(
        np.arange(0, sliceband),
        n_slices // sliceband
    ).reshape(1, 1, -1, 1)
    slice_times = tis_array + (slicedt * slice_numbers)
    # load images
    asl_img = Image(str(asl_name))
    t1_img = Image(str(t1_name))
    # check dimensions of t1 image to see if time series or not
    if t1_img.ndim == 3:
        t1_data = t1_img.data[..., np.newaxis]
    elif t1_img.ndim == 4:
        t1_data = t1_img.data
    # multiply asl sequence by satrecov model evaluated at TI
    numexp = - tis_array / t1_data
    num = 1 - np.exp(numexp)
    # divide asl sequence by satrecov model evaluated at actual slice-time
    denexp = - slice_times / t1_data
    den = 1 - np.exp(denexp)
    # evaluate scaling factors
    stcorr_factors = num / den
    stcorr_factors_img = Image(stcorr_factors, header=asl_img.header)
    # correct asl series
    stcorr_data = asl_img.data * stcorr_factors
    stcorr_img = Image(stcorr_data, header=asl_img.header)
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
    transforms = sorted(transform_dir.glob('**/*'))
    out_names = []
    for n, transform in enumerate(transforms):
        # naming intermediate file
        if n < 10:
            out_n = param_reg_name.parent / f'{param_reg_name.stem.split(".")[0]}_000{n}.nii.gz'
        else:
            out_n = param_reg_name.parent / f'{param_reg_name.stem.split(".")[0]}_00{n}.nii.gz'
        out_names.append(out_n)
        applyxfm(str(param_name), str(reffile), str(transform), str(out_n))
        cmd = [
            "applyxfm4D",
            str(param_name),
            str(reffile),
            str(out_n),
            str(transform),
            "-singlematrix"
        ]
        subprocess.run(cmd)
    # merge registered parameter volumes into one time series
    cmd = [
        'fslmerge',
        '-t',
        str(param_reg_name),
        *out_names
    ]
    subprocess.run(cmd)
    # remove intermediate file names
    for out_n in out_names:
        out_n.unlink()

def hcp_asl_moco(subject_dir, mt_factors, superfactor=1, cores=mp.cpu_count(), 
                 interpolation=3, nobandingcorr=False, outdir="hcp_asl"):
    """
    Full ASL correction and motion estimation pipeline.

    This function performs the full motion-correction pipeline for 
    the HCP ASL data. The steps of the pipeline are:
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
    #. Apply motion estimates to resulting T1 map so that it is 
    aligned with each of the volumes in the ASL series;
    #. Refined slice-timing correction

    Parameters
    ----------
    subject_dir : pathlib.Path
        Path to the subject's base directory.
    mt_factors : pathlib.Path
        Path to the pre-calculated MT correction scaling factors.
    superfactor : int, optional
        superfactor to use when using regtricks. Default is 1.
    cores : int, optional
        Number of cores regtricks will use. Default is the number 
        of cores available.
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
    """
    assert (isinstance(cores, int) and cores>0 and cores<=mp.cpu_count()), f"Number of cores should be an integer from 1-{mp.cpu_count()}."
    assert (isinstance(interpolation, int) and interpolation>=0 and interpolation<=5), "Order of interpolation should be an integer from 0-5."
    # asl sequence parameters
    ntis = 5
    iaf = "tc"
    ibf = "tis"
    tis = [1.7, 2.2, 2.7, 3.2, 3.7]
    rpts = [6, 6, 6, 10, 15]
    slicedt = 0.059
    sliceband = 10
    n_slices = 60
    # load json containing important file info
    json_dict = load_json(subject_dir/outdir)
    # original ASL series and bias field names
    asl_name = Path(json_dict['ASL_seq'])
    bias_name = json_dict['calib0_bias']
    calib_name = Path(json_dict['calib0_corr'])
    # create directories for results
    tis_dir_name = Path(json_dict['TIs_dir'])
    bcorr_dir = tis_dir_name / 'BiasCorr'
    distcorr_dir = tis_dir_name / 'DistCorr/FirstPass'
    moco_dir = tis_dir_name / 'MoCo'
    asln2m0_name = moco_dir / 'asln2m0.mat'
    satrecov_dir = tis_dir_name / 'SatRecov'
    create_dirs([
        tis_dir_name,
        bcorr_dir,
        distcorr_dir,
        moco_dir,
        asln2m0_name,
        satrecov_dir
    ])
    if not nobandingcorr:
        mtcorr_dir = tis_dir_name / 'MTCorr'
        stcorr_dir = tis_dir_name / 'STCorr'
        create_dirs([mtcorr_dir, stcorr_dir])

    # apply distortion corrections to the ASL series
    print("Applying distortion corrections to original ASL series.")
    gdc_name = Path(json_dict["ASL_dir"])/"gradient_unwarp/fullWarp_abs.nii.gz"
    gdc_warp = rt.NonLinearRegistration.from_fnirt(coefficients=str(gdc_name), 
                                                   src=str(asl_name), 
                                                   ref=str(asl_name), 
                                                   intensity_correct=True,
                                                   constrain_jac=(0.01, 100))
    dc_name = Path(json_dict["ASL_dir"])/"topup/WarpField_01.nii.gz"
    dc_warp = rt.NonLinearRegistration.from_fnirt(coefficients=str(dc_name),
                                                  src=str(asl_name),
                                                  ref=str(asl_name),
                                                  intensity_correct=True,
                                                  constrain_jac=(0.01, 100))
    gdc_dc_warp = rt.chain(gdc_warp, dc_warp)
    asl_gdc_dc = gdc_dc_warp.apply_to_image(src=str(asl_name),
                                            ref=str(asl_name),
                                            order=interpolation,
                                            cores=cores)
    asl_gdc_dc_name = distcorr_dir/"tis_gdc_dc.nii.gz"
    nb.save(asl_gdc_dc, asl_gdc_dc_name)

    # bias correct the ASL series
    print("Bias-correcting the distortion-corrected ASL series.")
    bcorr_img = bcorr_dir / 'tis_biascorr.nii.gz'
    fslmaths(str(asl_gdc_dc_name)).div(str(bias_name)).run(str(bcorr_img))

    # apply MT scaling factors to the bias-corrected ASL series
    if not nobandingcorr:
        print("MT-correcting the distortion-corrected ASL series.")
        mtcorr_name = mtcorr_dir / 'tis_mtcorr.nii.gz'
        mt_sfs = np.loadtxt(mt_factors)
        mt_arr = np.repeat(np.tile(mt_sfs, (86, 86, 1))[..., np.newaxis], 86, axis=-1)
        mt_img = nb.nifti1.Nifti1Image(mt_arr, affine=asl_gdc_dc.affine)
        biascorr_img = Image(str(bcorr_img))
        assert (len(mt_sfs) == biascorr_img.shape[2])
        mtcorr_img = Image(biascorr_img.data*mt_img.get_fdata(), header=biascorr_img.header)
        mtcorr_img.save(str(mtcorr_name))
        asl_corr = mtcorr_name
    else:
        asl_corr = bcorr_img

    # estimate satrecov model on distortion-, bias- and MT- corrected ASL series
    print("First satrecov model fit.")
    t1_name = _saturation_recovery(asl_corr, satrecov_dir, ntis, iaf, ibf, tis, rpts)
    t1_filt_name = _fslmaths_med_filter_wrapper(t1_name)

    # perform slice-time correction using estimated tissue params
    if not nobandingcorr:
        print("Performing initial ST correction.")
        stcorr_img, stfactors_img = _slicetiming_correction(mtcorr_name, t1_filt_name, tis, rpts, slicedt, sliceband, n_slices)
        stcorr_name = stcorr_dir / 'tis_stcorr.nii.gz'
        stcorr_img.save(stcorr_name)
        stfactors_name = stcorr_dir / 'st_scaling_factors.nii.gz'
        stfactors_img.save(stfactors_name)
        asl_corr = stcorr_name

    # register ASL series to calibration image
    print("Running mcflirt on calibration image and ASL series.")
    reg_name = moco_dir / 'initial_registration_TIs.nii.gz'
    mcflirt(str(asl_corr), reffile=json_dict['calib0_corr'], mats=True, out=str(reg_name))
    # rename mcflirt matrices directory
    orig_mcflirt = moco_dir / 'initial_registration_TIs.nii.gz.mat'
    if asln2m0_name.exists():
        shutil.rmtree(asln2m0_name)
    orig_mcflirt.rename(asln2m0_name)

    # register pre-ST-correction ASLn to ASL0
    print("Apply distortion and motion correction to original ASL series.")
    asln2m0_moco = rt.MotionCorrection.from_mcflirt(mats=str(asln2m0_name),
                                                    src=str(asl_name),
                                                    ref=json_dict['calib0_corr'])
    asln2asl0 = rt.chain(asln2m0_moco, asln2m0_moco.transforms[0].inverse())
    gdc_dc_asln2asl0 = rt.chain(gdc_dc_warp, asln2asl0)
    reg_gdc_dc = gdc_dc_asln2asl0.apply_to_image(str(asl_name), 
                                                 json_dict['calib0_corr'],
                                                 cores=cores,
                                                 order=interpolation)
    
    # apply moco to MT scaling factors image
    if not nobandingcorr:
        print("Apply motion estimates to the MT scaling factors image")
        mt_reg_img = asln2asl0.apply_to_image(src=mt_img, ref=mt_img, order=interpolation)
        print(mt_img.shape, mt_reg_img.shape)

    # apply bias-correction to motion- and distortion-corrected ASL series
    print("Apply bias correction to the distortion- and motion-corrected ASL series.")
    bias_img = nb.load(bias_name)
    reg_gdc_dc_biascorr = nb.nifti1.Nifti1Image(reg_gdc_dc.get_fdata()/bias_img.get_fdata()[..., np.newaxis],
                                                affine=reg_gdc_dc.affine)
    reg_gdc_dc_biascorr_name = moco_dir / 'reg_gdc_dc_tis_biascorr.nii.gz'
    nb.save(reg_gdc_dc_biascorr, reg_gdc_dc_biascorr_name)

    # apply MT-correction
    if not nobandingcorr:
        temp_reg_gdc_dc_mtcorr = moco_dir / 'temp_reg_dc_tis_mtcorr.nii.gz'
        reg_gdc_dc_mtcorr = nb.nifti1.Nifti1Image(reg_gdc_dc_biascorr.get_fdata()*mt_reg_img.get_fdata(),
                                                  affine=reg_gdc_dc.affine)
        nb.save(reg_gdc_dc_mtcorr, temp_reg_gdc_dc_mtcorr)
        asl_corr = temp_reg_gdc_dc_mtcorr
    else:
        asl_corr = reg_gdc_dc_biascorr_name

    # re-estimate satrecov model on distortion- and motion-corrected data
    print("Re-fitting the satrecov model since data has been motion-corrected.")
    satrecov_dir = tis_dir_name / 'SatRecov2'
    stcorr_dir = tis_dir_name / 'STCorr2'
    create_dirs([satrecov_dir, stcorr_dir])
    t1_name = _saturation_recovery(asl_corr, satrecov_dir, ntis, iaf, ibf, tis, rpts)
    t1_filt_name = _fslmaths_med_filter_wrapper(t1_name)
    
    if not nobandingcorr:
        # apply refined slice-time correction to registered- and distortion-corrected ASL series
        print("ST correcting the ASL series.")
        stcorr_img, stfactors_img = _slicetiming_correction(temp_reg_gdc_dc_mtcorr, 
                                                            t1_filt_name, 
                                                            tis, 
                                                            rpts, 
                                                            slicedt, 
                                                            sliceband, 
                                                            n_slices)
        stcorr_name = stcorr_dir / 'tis_stcorr.nii.gz'
        stfactors_name = stcorr_dir / 'st_scaling_factors.nii.gz'
        stcorr_img.save(str(stcorr_name))
        stfactors_img.save(str(stfactors_name))

        # combined MT and ST scaling factors
        print("Combining the ST and MT scaling factors into one set of scaling factors.")
        combined_factors_name = stcorr_dir / 'combined_scaling_factors.nii.gz'
        combined_factors_img = Image(stfactors_img.data*mt_reg_img.get_fdata(), header=stfactors_img.header)
        combined_factors_img.save(str(combined_factors_name))

        asl_corr = stcorr_name
    else:
        combined_factors_img = nb.nifti1.Nifti1Image(np.ones_like(reg_gdc_dc_biascorr.get_fdata()),
                                                     affine=reg_gdc_dc_biascorr.affine)
        combined_factors_name = moco_dir / 'combined_scaling_factors.nii.gz'
        nb.save(combined_factors_img, combined_factors_name)
    
    # save locations of important files in the json
    important_names = {
        'ASL_corr': str(asl_corr),
        'scaling_factors' : str(combined_factors_name)
    }
    update_json(important_names, json_dict)
    