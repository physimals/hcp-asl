from fsl.data import atlases
from fsl.data.image import Image
from fsl.wrappers import fslroi, fslmaths, LOAD
from fsl.wrappers.fnirt import invwarp, applywarp

import nibabel as nb

import numpy as np

import subprocess

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
    calib_dirs = [subject_dir/f'ASL/Calib/{c}' for c in ('Calib0', 'Calib1')]
    create_dirs(calib_dirs)

    # _A and _B directories
    a_dir = (subject_dir/f'{subid}_V1_A/scans').resolve(strict=True)
    b_dir = (subject_dir/f'{subid}_V1_B/scans').resolve(strict=True)

    # find T1 structural image
    t1_dir = list((a_dir).glob('**/*T1w'))[0]
    t1_dir = t1_dir/'resources/NIFTI/files'
    t1_name = list(t1_dir.glob(f'**/{subid}_*.nii.gz'))[0].resolve(strict=True)

    # create directories
    calib_dirs = [subject_dir/f'ASL/Calib/{c}' for c in ('Calib0', 'Calib1')]
    create_dirs(calib_dirs)
    aslt1_dir = t1_dir/"ASL"
    aslt1_dir.mkdir(exist_ok=True)

    # find mbpcasl sequence and sefms
    b_contents = b_dir.glob('**/*')
    for name in b_contents:
        if "mbPCASLhr_PA.nii.gz" in str(name):
            mbpcasl = name
        elif "PCASLhr_SpinEchoFieldMap_PA.nii.gz" in str(name):
            pa_sefm = name
        elif "PCASLhr_SpinEchoFieldMap_AP.nii.gz" in str(name):
            ap_sefm = name
            
    # split mbpcasl sequence into its calibration images
    calib_names = [d/f'calib{n}.nii.gz' for n, d in enumerate(calib_dirs)]
    [fslroi(str(mbpcasl), calib_name, n, 1) for calib_name, n in zip(calib_names, (88, 89))]

    # return helpful dictionary
    names = {
        "calib0_dir": calib_dirs[0],
        "calib0_name": calib_names[0],
        "calib1_dir": calib_dirs[1],
        "calib1_name": calib_names[1],
        "t1_dir": t1_dir,
        "t1_name": t1_name,
        "aslt1_dir": aslt1_dir,
        "pa_sefm": pa_sefm,
        "ap_sefm": ap_sefm
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
    seg = Image(np.where(pve.data>threshold, 1., 0.), header=pve.header)
    return seg

def generate_wmmask(aparc_aseg):
    """ 
    Generate binary WM mask in space of T1 image using FS aparc+aseg

    Args: 
        aparc_aseg: path to aparc_aseg in T1 space (not FS 256 1mm space!)

    Returns: 
        np.array logical WM mask in space of T1 image 
    """

    aseg_array = nb.load(aparc_aseg).get_data()
    wm = np.logical_or(aseg_array == 41, aseg_array == 2)
    return wm 

def linear_asl_reg(calib_name, results_dir, t1_name, t1_brain_name, wmpve_name):
    # get white matter segmentation for registration
    wmseg = binarise(wmpve_name)
    wmseg_name = results_dir/"wm_seg.nii.gz"
    wmseg.save(str(wmseg_name))
    # run asl_reg
    aslreg_cmd = [
        "asl_reg", 
        "-i", calib_name, "-o", results_dir, 
        "-s", t1_name, "--sbet", t1_brain_name, 
        "--tissseg", wmseg_name
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
    harvox = atlases.loadAtlas('harvardoxford-subcortical', 
                                         resolution=2.0)
    vent_img = Image(
        harvox.data[:,:,:,2] + harvox.data[:,:,:,13], 
        header=harvox.header
    )
    vent_img = fslmaths(vent_img).thr(0.1).bin().ero().run(LOAD)

    # apply inv(T1->MNI) registration
    t1_brain = (fslanatdir/"T1_biascorr_brain.nii.gz").resolve(strict=True)
    struct2mni = (fslanatdir/"T1_to_MNI_nonlin_coeff.nii.gz").resolve(strict=True)
    mni2struct = invwarp(str(struct2mni), str(t1_brain), LOAD)['out']
    vent_t1_img = applywarp(vent_img, str(t1_brain), LOAD, warp=mni2struct)['out']
    vent_t1_img = fslmaths(vent_t1_img).thr(0.9).bin().run(LOAD)
    return vent_t1_img
