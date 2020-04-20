"""
This contains a range of functions required to correct for the 
Magnetisation Transfer effect visible in the HCP data.
"""

import json
from pathlib import Path
from fsl.wrappers import fslmaths, LOAD, bet, fast
from initial_bookkeeping import create_dirs
import subprocess

def load_json(subject_dir):
    """
    Load json but with some error-checking to make sure it exists.
    If it doesn't exist, instruct user to run the first part of 
    the pipeline first.
    """
    json_name = subject_dir / 'ASL/ASL.json'
    if json_name.exists():
        with open(json_name, 'r') as infile:
            json_dict = json.load(infile)
    else:
        raise Exception(f'File {json_name} does not exist.' + 
                        'Please run initial_processing() first.')
    return json_dict

def correct_M0(subject_dir, mt_factors):
    """
    Correct the M0 images for a particular subject whose data 
    is stored in `subject_dir`. The corrections to be 
    performed include:
        - Bias-field correction
        - Magnetisation Transfer correction
    
    Inputs
        - `subject_dir` = pathlib.Path object specifying the 
            subject's base directory
        - `mt_factors` = pathlib.Path object specifying the 
            location of empirically estimated MT correction 
            scaling factors
    """
    # load json containing info on where files are stored
        # if doesn't exist, throw error and ask user to run
        # initial_processing() from initial_bookkeeping.py first
    load_json(subject_dir)
    
    # do for both m0 images for the subject, calib0 and calib1
    calib_names = [json_dict['CALIB0_img'], json_dict['CALIB1_img']]
    for calib_name in calib_names:
        # run BET on m0 image
        betted_m0 = bet(calib_name, LOAD)
        # estimate bias field on brain-extracted m0 image
            # create FAST directory
        fast_dir = Path(calib_name).parent / 'FAST'
        create_dirs(fast_dir)
            # run FAST, storing results in directory
        fast_base = fast_dir / 'FAST'
        fast(
            betted_m0['output'], # output of bet
            out=LOAD, 
            t=3, # image type, 3=PD image
            B=True, # output bias-corrected image
            b=True # output estimated bias field
        )
        # apply bias field to original m0 image (i.e. not BETted)

        # apply mt_factors to bias-corrected m0 image

        # save bias field, bias-corrected m0 and mt-corrected m0

        # add locations of above files to the json
    pass