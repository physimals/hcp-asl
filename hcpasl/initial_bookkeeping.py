"""
A collection of functions to be used for the inital set-up for a 
subject. These set-up tasks include:

    - Creating ASL sub-directories in the subject's directory
    - Finding T1w, fieldmap and mbPCASL directories and scans
    - Splitting mbPCASL sequence into its components of tis, 
        and calibration images
    - Creating a json to keep track of import files and 
        directories
"""

from pathlib import Path
from fsl.wrappers.misc import fslroi
from fsl.wrappers.fsl_anat import fsl_anat
import json

def create_dirs(dir_list, parents=True, exist_ok=True):
    """
    Given a list of pathlib.Path directories, `dir_list`, create 
    these directories.

    Default behaviour is to create parent directories if these 
    don't yet exist and to not throw an error if a directory 
    already exists.

    Inputs:
        - `dir_list` = list of pathlib.Path objects which are 
            the directories you wish to create
        - `parents` = Path.mkdir() argument which, if `True`, will 
            create parent directories if they do not yet exist. 
            See pathlib.Path.mkdir() for more details. Default 
            here is `True`.
        - `exist_ok` = Path.mkdir() argument which, if `True`, 
            doesn't throw an error if the directory already exists.
            See pathlib.Path.mkdir() for more details. Default 
            here is `True`.
    """
    for directory in dir_list:
        directory.mkdir(parents=parents, exist_ok=exist_ok)

def initial_processing(subject_dir, mbpcasl, structural, surfaces):
    """
    Perform initial processing for the subject directory provided.
    These initial processing includes:

    - Creating ASL sub-directories in the subject's directory
    - Finding T1w, fieldmap and mbPCASL directories and scans
    - Splitting mbPCASL sequence into its components of tis, 
        and calibration images
    - Run fsl_anat on the subject's structural image
    - Creating a json to keep track of import files and 
        directories
    
    Input:
        - `subject_dir` = a pathlib.Path object for the subject's
            data directory
    """
    # get subject name
    subject_name = subject_dir.parts[-1]

    # create ${subject_dir}/ASL and ${subject_dir}/T1w/Results/ASL 
    # directories
    asl_dir = subject_dir / 'ASL'
    tis_dir = asl_dir / 'TIs'
    calib_dir = asl_dir / 'Calib'
    calib0_dir = calib_dir / 'Calib0'
    calib1_dir = calib_dir / 'Calib1'
    strucasl_dir = subject_dir / 'T1w/ASL'
    create_dirs([asl_dir, tis_dir, calib0_dir, calib1_dir, strucasl_dir])

    # find sub-directories
    # structural
    t1_dir = subject_dir / 'T1w'
    t1_name = structural['struct']
    t1_brain_name = structural['sbrain']
    
    # output names
    tis_name = tis_dir / 'tis.nii.gz'
    calib0_name = calib0_dir / 'calib0.nii.gz'
    calib1_name = calib1_dir / 'calib1.nii.gz'
    # get tis
    fslroi(str(mbpcasl), tis_name, 0, 86)
    # get calibration images
    fslroi(str(mbpcasl), calib0_name, 88, 1)
    fslroi(str(mbpcasl), calib1_name, 89, 1)

    # add filenames to a dictionary to be saved to a json
    json_name = asl_dir / 'ASL.json'
    fields = [
        "T1w_dir",
        "T1w_acpc",
        "T1w_acpc_brain",
        "ASL_seq",
        "ASL_dir",
        "TIs_dir",
        "structasl",
        "calib_dir",
        "calib0_dir",
        "calib1_dir",
        "calib0_img",
        "calib1_img",
        "L_mid",
        "R_mid",
        "L_pial",
        "R_pial",
        "L_white",
        "R_white",
        "json_name"
    ]
    field_values = [
        t1_dir,
        t1_name,
        t1_brain_name,
        tis_name,
        asl_dir,
        tis_dir,
        strucasl_dir,
        calib_dir,
        calib0_dir,
        calib1_dir,
        calib0_name,
        calib1_name,
        surfaces['L_mid'],
        surfaces['R_mid'],
        surfaces['L_pial'],
        surfaces['R_pial'],
        surfaces['L_white'],
        surfaces['R_white'],
        json_name
    ]
    names_dict = {}
    for key, value in zip(fields, field_values):
        names_dict[key] = str(value)
    with open(json_name, 'w') as fp:
        json.dump(names_dict, fp, sort_keys=True, indent=4)