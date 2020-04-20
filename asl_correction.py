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

    # perform initial slice-timing correction using estimated tissue params
        # median filter the parameter estimates

    # motion estimation from ASL to M0 image

    # apply inverse transformations to parameter estimates to align them with 
    # the individual frames of the ASL series

    # second slice-timing correction using registered parameter estimates

    # apply motion estimates to slice-timing corrected ASL series

    # save locations of important files in the json
    pass