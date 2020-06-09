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

from .initial_bookkeeping import create_dirs
from .m0_mt_correction import load_json, update_json
from fsl.wrappers import fslmaths, LOAD
from fsl.wrappers.flirt import mcflirt, applyxfm, applyxfm4D
from fsl.data.image import Image
from fabber import Fabber, percent_progress
import sys
from pathlib import Path
import shutil
import subprocess
import numpy as np
import regtricks as rt
def _satrecov_worker(control_name, satrecov_dir, tis, rpts, ibf, spatial):
    """
    Runs fabber's saturation recovery model on the given sequence 
    of control images.

    Inputs:
        - `control_name` = name of control sequence
        - `satrecov_dir` = parent directory for the satrecov 
            results. Results from this will be stored either 
            in {`satrecov_dir`}/spatial or 
            {`satrecov_dir`}/nospatial depending on value of 
            `spatial`
        - `tis` = list of TIs for the ASL sequence
        - `rpts` = list of repeats for each TI in the sequence
        - `ibf` = input format of the sequence
        - `spatial` = Boolean for whether to run fabber in 
            spatial (True) or non-spatial (False) mode
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
    Given and ASL time series, `asl_name`, and the sequence details, 
    split the series into 2 files, even and odd. Save these with the 
    base of `asl_name`_even.nii.gz for example.

    Inputs:
        - `asl_name` = pathlib.Path object for the ASL series to be 
            split
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
    Wrapper function for Fabber's `satrecov` model.

    Given an ASL time series, `asl_name`, estimate Fabber's 
    `satrecov` model on the series' control images. This 
    function will assume the ASL data is HCP-ASL data which 
    is ###insert tag-control / control-tag here###.

    Fabber will be called twice:
        - once with spatial mode off
        - once with spatial mode on, initialised from the 
            prior run
    
    Inputs:
        - `asl_name` = pathlib.Path object for the ASL series on 
            which the model will be estimated
        - `results_dir` = pathlib.Path object in which to store the 
            `nospatial` and `spatial` parameter estimates
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
    Performs rescaling of ASL series, `asl_name`, accounting for 
    the estimated T1t in a given voxel, `t1_name`.

    satrecov_model: S(t) = M_0 * (1 - exp{-t/T1t})

    Divides the value in a voxel by the above model evaluated at 
    t = TI + n*slicedt where n is the slice number of the voxel 
    within a band. This accounts for the voxel being imaged at a 
    different time than its supposed TI.

    Multiplies the value in a voxel by the above model evaluated 
    at t = TI, i.e. scales the values as if they had been imaged 
    at the TI that was specified in the ASL sequence.

    Returns the slice-timing corrected ASL series, `stcorr_img`, 
    and the scaling factors used to perform the correction, 
    `stcorr_factors_img` where:
        `stcorr_factors` = S(TI) / S(TI + n*slicedt)
        `stcorr_img` = `stcorr_factors` * `asl_name`
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
    Given a parameter map, `param_name`, and a series of motion 
    estimates, `transform_dir`, apply the motion estimates to 
    the parameter map and obtain a time series of the map 
    in the frame described by the motion estimates.

    Inputs:
        - `param_name` = pathlib.Path object for the parameter 
            estimate
        - `transform_dir` = pathlib.Path object for the 
            directory containing mcflirt motion estimates
        - `reffile` = filename for reference image in mcflirt 
            motion estimate
        - `param_reg_name` = name of output file
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

def hcp_asl_moco(subject_dir, mt_factors, superlevel=1, cores=1, order=3):
    """
    This function performs the full motion-correction pipeline for 
    the HCP ASL data. The steps of the pipeline include:
    - Bias-field correction
    - MT correction
    - Saturation recovery
    - Initial slice-timing correction
    - Motion estimation
    - Second slice-timing correction
    - Registration

    Inputs
        - `subject_dir` = pathlib.Path object specifying the 
            subject's base directory
        - `mt_factors` = pathlib.Path object specifying the 
            location of empirically estimated MT correction 
            scaling factors
    """
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
    json_dict = load_json(subject_dir)
    # create directories for results
    tis_dir_name = Path(json_dict['TIs_dir'])
    first_pass_dir = tis_dir_name / 'FirstPass'
    second_pass_dir = tis_dir_name / 'SecondPass'
    create_dirs([tis_dir_name, first_pass_dir, second_pass_dir])
    # original ASL series and bias field names
    asl_name = Path(json_dict['ASL_seq'])
    bias_name = json_dict['calib0_bias']
    old_m02asl = first_pass_dir / 'MoCo/m02asln.mat'
    # iterate over first and second passes
    for n, iteration in enumerate((first_pass_dir, second_pass_dir)):
        bcorr_dir = iteration / 'BiasCorr'
        mtcorr_dir = iteration / 'MTCorr'
        satrecov_dir = iteration / 'SatRecov'
        stcorr_dir = iteration / 'STCorr'
        moco_dir = iteration / 'MoCo'
        asln2m0_name = moco_dir / 'asln2m0.mat'
        m02asln_name = moco_dir / 'm02asln.mat'
        asln2asl0_name = moco_dir / 'asln2asl0.mat'
        asl02asln_name = moco_dir / 'asl02asln.mat'
        create_dirs([
            bcorr_dir,
            mtcorr_dir,
            satrecov_dir,
            stcorr_dir,
            moco_dir,
            asln2m0_name,
            m02asln_name,
            asln2asl0_name,
            asl02asln_name
        ])
        # bias correct the original ASL series
        bcorr_img = bcorr_dir / 'tis_biascorr.nii.gz'
        if n == 1:
            # register bias field to ASL series
            reg_bias_name = bcorr_dir / 'bias_reg.nii.gz'
            _register_param(bias_name, old_m02asl, bias_name, reg_bias_name)
            old_m02asl = rt.MotionCorrection(old_m02asl)
            Image(
                old_m02asl.apply_to_image(
                    bias_name,
                    bias_name,
                    superlevel=superlevel,
                    cores=cores,
                    order=order
                )
            ).save(str(reg_bias_name))
            bias_name = reg_bias_name
        fslmaths(str(asl_name)).div(str(bias_name)).run(str(bcorr_img))
        # apply MT scaling factors to the bias-corrected ASL series
        mtcorr_name = mtcorr_dir / 'tis_mtcorr.nii.gz'
        fslmaths(str(bcorr_img)).mul(str(mt_factors)).run(str(mtcorr_name))
        # estimate satrecov model on bias and MT corrected ASL series
        t1_name = _saturation_recovery(mtcorr_name, satrecov_dir, ntis, iaf, ibf, tis, rpts)
        t1_filt_name = _fslmaths_med_filter_wrapper(t1_name)
        # perform slice-time correction using estimated tissue params
        stcorr_img, stfactors_img = _slicetiming_correction(mtcorr_name, t1_filt_name, tis, rpts, slicedt, sliceband, n_slices)
        stcorr_name = stcorr_dir / 'tis_stcorr.nii.gz'
        stcorr_img.save(stcorr_name)
        stfactors_name = stcorr_dir / 'st_scaling_factors.nii.gz'
        stfactors_img.save(stfactors_name)
        # register ASL series to calibration image
        reg_name = moco_dir / 'initial_registration_TIs.nii.gz'
        mcflirt(stcorr_img, reffile=json_dict['calib0_mc'], mats=True, out=str(reg_name))
        # rename mcflirt matrices directory
        orig_mcflirt = moco_dir / 'initial_registration_TIs.nii.gz.mat'
        if asln2m0_name.exists():
            shutil.rmtree(asln2m0_name)
        orig_mcflirt.rename(asln2m0_name)
        # get motion estimates from ASLn to ASL0 (and their inverses)
        asl2m0_list = sorted(asln2m0_name.glob('**/MAT*'))
        m02asl0 = np.linalg.inv(np.loadtxt(asl2m0_list[0]))
        for n, xform in enumerate(asl2m0_list):
            if n == 0:
                fwd_xform = np.eye(4)
            else:
                fwd_xform = m02asl0 @ np.loadtxt(xform)
            inv_xform = np.linalg.inv(fwd_xform)
            np.savetxt(m02asln_name / xform.stem, np.linalg.inv(np.loadtxt(xform)))
            np.savetxt(asln2asl0_name / xform.stem, fwd_xform)
            np.savetxt(asl02asln_name / xform.stem, inv_xform)
    # register pre-ST-correction ASLn to ASL0
    temp_reg_mtcorr = moco_dir / 'temp_reg_tis_mtcorr.nii.gz'
    asln2m0_moco = rt.MotionCorrection(asln2m0_name)
    asln2asl0 = rt.chain(asln2m0_moco, asln2m0_moco.transforms[0].inverse())
    reg_mtcorr = Image(asln2asl0.apply_to_image(
        str(mtcorr_name), 
        json_dict['calib0_mc'],
        superlevel=superlevel,
        cores=cores,
        order=order
    ))
    reg_mtcorr.save(str(temp_reg_mtcorr))

    # estimate satrecov model on motion-corrected data
    satrecov_dir = iteration / 'SatRecov2'
    stcorr_dir = iteration / 'STCorr2'
    create_dirs([satrecov_dir, stcorr_dir])
    t1_name = _saturation_recovery(temp_reg_mtcorr, satrecov_dir, ntis, iaf, ibf, tis, rpts)
    t1_filt_name = _fslmaths_med_filter_wrapper(t1_name)
    # apply asl0 to asln registrations to new t1 map
    reg_t1_filt_name = t1_filt_name.parent / f'{t1_filt_name.stem.split(".")[0]}_reg.nii.gz'
    reg_t1_filt = Image(
        asln2asl0.inverse().apply_to_image(
            str(t1_filt_name),
            json_dict['calib0_mc'],
            superlevel=superlevel,
            cores=cores,
            order=order
        )
    )
    reg_t1_filt.save(str(reg_t1_filt_name))
    # perform slice-time correction using estimated tissue params
    stcorr_img, stfactors_img = _slicetiming_correction(mtcorr_name, reg_t1_filt_name, tis, rpts, slicedt, sliceband, n_slices)
    # save images
    stcorr_name = stcorr_dir / 'tis_stcorr.nii.gz'
    stfactors_name = stcorr_dir / 'st_scaling_factors.nii.gz'
    stcorr_img.save(str(stcorr_name))
    stfactors_img.save(str(stfactors_name))
    # combined MT and ST scaling factors
    combined_factors_name = stcorr_dir / 'combined_scaling_factors.nii.gz'
    fslmaths(str(stfactors_name)).mul(str(mt_factors)).run(str(combined_factors_name))
    # save locations of important files in the json
    important_names = {
        'ASL_stcorr': str(stcorr_name),
        'scaling_factors': str(combined_factors_name)
    }
    update_json(important_names, json_dict)