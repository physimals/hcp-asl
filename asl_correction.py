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

from m0_mt_correction import load_json
from fsl.wrappers import fslmaths, LOAD
from fsl.wrappers.flirt import mcflirt, applyxfm
from fsl.data.image import Image
from fabber import Fabber, percent_progress
from initial_bookkeeping import create_dirs
from pathlib import Path
import subprocess
def _satrecov_worker():
    pass

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

def _saturation_recovery(asl_name, results_dir):
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
    # asl sequence parameters
    ntis = 5
    iaf = "tc"
    ibf = "tis"
    tis = [1.7, 2.2, 2.7, 3.2, 3.7]
    rpts = [6, 6, 6, 10, 15]
    # obtain control images of ASL series
    _split_tag_control(asl_name, ntis, iaf, ibf, rpts)
    # satrecov nospatial
    # satrecov spatial
    # return T1t map?
    pass

def hcp_asl_moco(subject_dir, mt_factors):
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
    """
    # load json containing important file info
    json_dict = load_json(subject_dir)

    # create directories for results
    asl_dir_name = Path(json_dict['ASL_dir'])
    biascorr_dir_name = asl_dir_name / 'BIASCORR'
    mtcorr_dir_name = asl_dir_name / 'MTCORR'
    satrecov_dir_name = asl_dir_name / 'SATRECOV'
    moco_dir_name = asl_dir_name / 'MOCO'
    # this is getting messy, maybe a function to make
    # the directories given a parent and a list of 
    # strings to create as sub-dirs within the parent?
    create_dirs([
        biascorr_dir_name, 
        mtcorr_dir_name, 
        satrecov_dir_name,
        moco_dir_name
    ])

    # bias-correction of original ASL series
    asl_name = json_dict['ASL_seq']
    bias_name = json_dict['calib0_bias'] # possibly take mean of both bias estimates?
    # possibly some rough registration from M0 to mean of ASL series
        # if doing the above, is it worth running BET on M0 images again
        # and changing f parameter so that the brain-mask is larger?
    biascorr_name = biascorr_dir_name / 'tis_biascorr.nii.gz'
    fslmaths(str(asl_name)).div(str(bias_name)).run(str(biascorr_name))
    
    # apply MT scaling factors to the bias-corrected ASL series
    mtcorr_name = mtcorr_dir_name / 'tis_mtcorr.nii.gz'
    fslmaths(str(biascorr_name)).mul(str(mt_factors)).run(str(mtcorr_name))

    # estimate satrecov model on bias-corrected, MT-corrected ASL series
    _saturation_recovery(mtcorr_name, satrecov_dir_name)
    # perform initial slice-timing correction using estimated tissue params
        # median filter the parameter estimates

    # motion estimation from ASL to M0 image

    # apply inverse transformations to parameter estimates to align them with 
    # the individual frames of the ASL series

    # second slice-timing correction using registered parameter estimates

    # apply motion estimates to slice-timing corrected ASL series

    # save locations of important files in the json