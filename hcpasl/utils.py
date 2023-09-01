import json
import logging
import shutil
import subprocess
from importlib.resources import path as resource_path
from pathlib import Path
import os

import numpy as np
from fsl.data import atlases
from fsl.data.image import Image
from fsl.wrappers import LOAD, fslmaths, fslroi
from fsl.wrappers.fnirt import applywarp, invwarp
from fsl.wrappers.misc import fslroi

from . import resources

ASL_SHAPE = (86, 86, 60, 86)


def create_dirs(dir_list, parents=True, exist_ok=True):
    """
    Creates directories in a list.

    Default behaviour is to create parent directories if these
    don't yet exist and not throw an error if a directory
    already exists.

    Parameters
    ----------
    dir_list : list of pathlib.Path objects
        The directories to be created.
    parents : bool
        Create parent directories if they do not yet exist.
        See pathlib.Path.mkdir() for more details. Default
        here is `True`.
    exist_ok : bool
        Don't throw an error if the directory already exists.
        See pathlib.Path.mkdir() for more details. Default
        here is `True`.
    """
    for directory in dir_list:
        directory.mkdir(parents=parents, exist_ok=exist_ok)


def load_json(subject_dir):
    """
    Load json but with some error-checking to make sure it exists.
    If it doesn't exist, instruct user to run the first part of
    the pipeline first.

    Parameters
    ----------
    subject_dir : pathlib.Path
        Path to subject's base directory.

    Returns
    -------
    dict
        Dictionary containing important file paths.
    """
    json_name = subject_dir / "ASL/ASL.json"
    if json_name.exists():
        with open(json_name, "r") as infile:
            json_dict = json.load(infile)
    else:
        raise Exception(
            f"File {json_name} does not exist."
            + "Please run initial_processing() first."
        )
    return json_dict


def update_json(new_dict, old_dict):
    """
    Add key-value pairs to the dictionary.

    Add the key-value pairs in `new_dict` to `old_dict`
    and save the resulting dictionary to the json found in
    `old_dict['json_name']`.

    Parameters
    ----------
    new_dict : dict
        New key-value pairs to add to the json of important
        file names.
    old_dict : dict
        Dictionary to be updated. Also has a field containing
        the location of the json to update.
    """
    old_dict.update(new_dict)
    with open(Path(old_dict["json_name"]), "w") as fp:
        json.dump(old_dict, fp, sort_keys=True, indent=4)


def setup(subject_dir):
    """
    Perform the initial set up for the MT Estimation pipeline.

    The setup includes finding necessary files (mbPCASL sequence,
    T1 structural image, spin echo field map images), creating
    directories and extracting the calibration images from the
    mbPCASL sequence.

    Parameters
    ----------
    subject_dir: pathlib.Path
        Path to the subject's directory

    Returns
    -------
    names: dict
        Dictionary of important filenames which will be used in the
        rest of the pipeline.
    """
    # make sure subject_dir exists
    subject_dir.resolve(strict=True)
    subid = subject_dir.parts[-1]

    # create directories
    calib_dirs = [subject_dir / f"ASL/Calib/{c}" for c in ("Calib0", "Calib1")]
    create_dirs(calib_dirs)

    # mbPCASLhr_unproc and Structural_preproc directories
    mbpcasl_dir = subject_dir / f"{subid}_V1_MR/resources/mbPCASLhr_unproc/files"
    struct_dir = subject_dir / f"{subid}_V1_MR/resources/Structural_preproc/files"

    # find T1 structural image
    t1_dir = struct_dir / f"{subid}_V1_MR/T1w"
    t1, t1_brain = [
        (t1_dir / f"T1w_acpc_dc_restore{suf}.nii.gz").resolve(strict=True)
        for suf in ("", "_brain")
    ]
    aparc_aseg, ribbon, wmparc = [
        (t1_dir / f"{name}.nii.gz").resolve(strict=True)
        for name in ("aparc+aseg", "ribbon", "wmparc")
    ]

    # create directories
    calib_dirs = [subject_dir / f"ASL/Calib/{c}" for c in ("Calib0", "Calib1")]
    create_dirs(calib_dirs)
    aslt1_dir = subject_dir / "T1wASL"
    aslt1_dir.mkdir(exist_ok=True, parents=True)

    # find mbpcasl sequence and sefms
    mbpcasl = (mbpcasl_dir / f"{subid}_V1_MR_mbPCASLhr_PA.nii.gz").resolve(strict=True)
    pa_sefm, ap_sefm = [
        (mbpcasl_dir / f"{subid}_V1_MR_PCASLhr_SpinEchoFieldMap_{suf}.nii.gz").resolve(
            strict=True
        )
        for suf in ("PA", "AP")
    ]

    # split mbpcasl sequence into its calibration images
    calib_names = [d / f"calib{n}.nii.gz" for n, d in enumerate(calib_dirs)]
    [
        fslroi(str(mbpcasl), calib_name, n, 1)
        for calib_name, n in zip(calib_names, (88, 89))
    ]

    # return helpful dictionary
    names = {
        "calib0_dir": calib_dirs[0],
        "calib0_name": calib_names[0],
        "calib1_dir": calib_dirs[1],
        "calib1_name": calib_names[1],
        "t1_dir": t1_dir,
        "t1_name": t1,
        "t1_brain_name": t1_brain,
        "aparc_aseg": aparc_aseg,
        "ribbon": ribbon,
        "wmparc": wmparc,
        "aslt1_dir": aslt1_dir,
        "pa_sefm": pa_sefm,
        "ap_sefm": ap_sefm,
    }
    return names


def binarise(pve_name, threshold=0.5):
    """
    Binarise a PVE image given a threshold.

    The default threshold is 0.5, the threshold used by asl_reg
    to obtain a white matter segmentation from a white matter
    PVE image.
    """
    pve = Image(str(pve_name))
    seg = Image(np.where(pve.data > threshold, 1.0, 0.0), header=pve.header)
    return seg


def linear_asl_reg(calib_name, results_dir, t1_name, t1_brain_name, wm_mask):
    # run asl_reg
    aslreg_cmd = [
        "asl_reg",
        "-i",
        calib_name,
        "-o",
        results_dir,
        "-s",
        t1_name,
        "--sbet",
        t1_brain_name,
        "--tissseg",
        wm_mask,
    ]
    subprocess.run(aslreg_cmd, check=True)


def get_ventricular_csf_mask(fslanatdir, interpolation=3):
    """
    Get a ventricular mask in T1 space.

    Register the ventricles mask from the harvardoxford-subcortical
    2mm atlas to T1 space, using the T1_to_MNI_nonlin_coeff.nii.gz
    in the provided fsl_anat directory.

    Parameters
    ----------
    fslanatdir: pathlib.Path
        Path to an fsl_anat output directory.
    interpolation: int
        Order of interpolation to use when performing registration.
        Default is 3.
    """
    # get ventricles mask from Harv-Ox
    atlases.rescanAtlases()
    harvox = atlases.loadAtlas("harvardoxford-subcortical", resolution=2.0)
    vent_img = Image(
        harvox.data[:, :, :, 2] + harvox.data[:, :, :, 13], header=harvox.header
    )
    vent_img = fslmaths(vent_img).thr(0.1).bin().ero().run(LOAD)

    # apply inv(T1->MNI) registration
    t1_brain = (fslanatdir / "T1_biascorr_brain.nii.gz").resolve(strict=True)
    struct2mni = (fslanatdir / "T1_to_MNI_nonlin_coeff.nii.gz").resolve(strict=True)
    mni2struct = invwarp(str(struct2mni), str(t1_brain), LOAD)["out"]
    vent_t1_img = applywarp(vent_img, str(t1_brain), LOAD, warp=mni2struct)["out"]
    vent_t1_img = fslmaths(vent_t1_img).thr(0.9).bin().run(LOAD)
    return vent_t1_img


def split_mbpcasl(mbpcasl, tis_name, calib0_name, calib1_name):
    """
    Split mbPCASLhr sequence into its constituent ASL series and
    calibration images.
    """
    fslroi(str(mbpcasl), str(tis_name), 0, 86)
    fslroi(str(mbpcasl), str(calib0_name), 88, 1)
    fslroi(str(mbpcasl), str(calib1_name), 89, 1)


def setup_logger(file_path):
    """
    Convenience function which returns a logger object with a
    specified name and level of reporting.

    Parameters
    ----------
    logger_name : str
        Name of the logger
    out_name : str
        Desired path for the logfile
    level : str (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        Desired level of reporting for the log. See Python's logging
        documentation for more info.
    verbose : bool, default=False
        If True, information will also be sent to the terminal via a
        StreamHandler as well as to a log file.
    mode : str, default="w"
        The mode of operation for the FileHandler. The default mode,
        "w", overwrites a logfile of the same name if it exists.
    """

    # set up logger's base reporting level and formatting
    logger = logging.getLogger()
    logger.setLevel("INFO")
    formatter = logging.Formatter(
        "%(asctime)s %(levelname)s %(module)s/%(funcName)s: %(message)s"
    )

    # set up FileHandler and StreamHanlder
    handlers = [logging.FileHandler(file_path, mode="w")]
    handlers.append(logging.StreamHandler())

    # add formatting to handlers and add to logger
    for handler in handlers:
        handler.setFormatter(formatter)
        logger.addHandler(handler)


def get_package_data_name(name):
    """
    This function returns the filepath to a given file which are
    included as part of the pipeline.

    The choices are:
    * report_template.ipynb
    * scaling_factors.txt
    """
    p = resource_path(resources, name)
    with p as filename:
        name = filename
    return name


def copy_oxford_asl_inputs(input_dict, output_dir):
    """
    Take a dictionary of inputs to oxford_asl and copy them into one directory.
    """
    # create new output directory
    output_dir.mkdir(exist_ok=True, parents=True)

    # create dict of new filenames
    new_dict = {key: output_dir / val.name for key, val in input_dict.items()}

    # copy files from original locations in input_dict to new location in new_dict
    for old_loc, new_loc in zip(input_dict.values(), new_dict.values()):
        shutil.copy2(old_loc, new_loc)


def make_motion_fov_mask(mc_transform, src, ref, cores=1):
    """Generate a FoV array in the reference of the motion transform
    Invert the transform and map the FoV into the source space
    Apply logical all across time to get valid FoV voxels
    Save mask in source space"""

    fov = np.ones(ref.size)
    fov_motion = mc_transform.apply_to_array(fov, ref, src, order=1, cores=cores)
    fov_valid = (fov_motion > 0.9).all(-1)
    return src.make_nifti(fov_valid)


def subprocess_popen(cmd, **kwargs):
    logging.info(cmd)
    with subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, **kwargs
    ) as process:
        while 1:
            retcode = process.poll()
            line = process.stdout.readline().decode("utf-8").replace("\n", "")
            if line:
                logging.info(line)
            if retcode is not None:
                break
        if retcode != 0:
            msg = f"Subprocess {cmd} failed with exit code {retcode}."
            logging.exception(msg)
            raise RuntimeError(msg)


def get_roi_stats_script():
    """FSL script can have two different names (with or without .py extension)"""
    roi_script_name = Path(os.environ["FSLDIR"]) / "bin/oxford_asl_roi_stats"
    if not roi_script_name.exists():
        roi_script_name = roi_script_name.with_suffix(".py")
    if not roi_script_name.exists():
        raise RuntimeError("Cannot find oxford_asl_roi_stats within $FSLDIR")
    return roi_script_name
