"""
Perform the setup necessary for estimating the MT Effect. This 
includes finding necessary files, creating results directories 
and running fsl_anat on the structural image.
"""

from pathlib import Path
from fsl.wrappers.misc import fslroi
from fsl.wrappers.fsl_anat import fsl_anat
from fsl.wrappers.flirt import applyxfm
from fsl.wrappers import fslmaths, LOAD, bet, fast
from fsl.data.image import Image
import json
from initial_bookkeeping import create_dirs
import subprocess

PVE_NAMES = {
    'csf': 'T1_fast_pve_0.nii.gz',
    'gm': 'T1_fast_pve_1.nii.gz',
    'wm': 'T1_fast_pve_2.nii.gz'
}

def setup_mtestimation(subject_dir, rois=['wm',]):
    """
    Perform the initial processing needed for estimation of the 
    MT Effect. This includes:
    - Creating sub-directories for storing the results
    - Finding T1 and mbPCASL directories and scans
    - Split mbPCASL sequence into its constituent components
    - Run fsl_anat
    - Create a json to keep track of important files and 
        directories
    """
    # get subject name
    subject_name = subject_dir.parts[-1]
    print(subject_name)

    # create results directories
    asl_dir = subject_dir / 'ASL'
    calib0_dir = asl_dir / 'Calib/Calib0'
    calib1_dir = asl_dir / 'Calib/Calib1'
    create_dirs([asl_dir, calib0_dir, calib1_dir])

    # obtain calibration images from ASL sequence
    mbpcasl_dir = list((subject_dir / f'{subject_name}_V1_B/scans').glob('**/*mbPCASLhr'))[0]
    mbpcasl = mbpcasl_dir / f'resources/NIFTI/files/{subject_name}_V1_B_mbPCASLhr_PA.nii.gz'
    calib0_name = calib0_dir / 'calib0.nii.gz'
    calib1_name = calib1_dir / 'calib1.nii.gz'
    fslroi(str(mbpcasl), calib0_name, 88, 1)
    fslroi(str(mbpcasl), calib1_name, 89, 1)

    # initialise dict
    json_name = asl_dir / 'ASL.json'
    important_dict = {
        "calib_dir": str(calib0_dir.parent),
        "calib0_dir": str(calib0_dir),
        "calib1_dir": str(calib1_dir),
        "calib0_img": str(calib0_name),
        "calib1_img": str(calib1_dir),
        "json_name": str(json_name)
    }

    # structural directory
    t1_dir = list((subject_dir / f'{subject_name}_V1_A/scans').glob('**/*T1w'))[0]
    struc_dir = t1_dir / 'resources/NIFTI/files'
    struc_name = list(struc_dir.glob(f'**/{subject_name}_*.nii.gz'))[0]
    fsl_anat_dir = struc_dir / 'ASL/struc'
    calib0struct_dir = struc_dir / 'ASL/Calib/Calib0'
    calib1struct_dir = struc_dir / 'ASL/Calib/Calib1'
    create_dirs([calib0struct_dir, calib1struct_dir])

    # run fsl_anat
    fsl_anat(str(struc_name), str(fsl_anat_dir), clobber=True, nosubcortseg=True)
    fsl_anat_dir = fsl_anat_dir.parent / f'{fsl_anat_dir.stem}.anat'
    t1_name = fsl_anat_dir / 'T1_biascorr.nii.gz'
    t1_brain_name = fsl_anat_dir / 'T1_biascorr_brain.nii.gz'
    # bias-field correction
    for calib_name in (calib0_name, calib1_name):
        calib_name_stem = calib_name.stem.split('.')[0]
        # run bet
        betted_m0 = bet(str(calib_name), LOAD)
        # create directories to save results
        fast_dir = calib_name.parent / 'FAST'
        biascorr_dir = calib_name.parent / 'BiasCorr'
        create_dirs([fast_dir, biascorr_dir])
        # run FAST on brain-extracted m0 image
        fast_base = fast_dir / calib_name_stem
        fast(
            betted_m0['output'],
            out=str(fast_base),
            type=3,
            b=True,
            nopve=True
        )
        bias_name = fast_dir / f'{calib_name_stem}_bias.nii.gz'
        # apply bias field to original m0 image (i.e. not BETted)
        biascorr_name = biascorr_dir / f'{calib_name_stem}_restore.nii.gz'
        fslmaths(str(calib_name)).div(str(bias_name)).run(str(biascorr_name))
        # obtain registration from structural to calibration image
        mask_dir = biascorr_name.parent / 'masks'
        create_dirs([mask_dir, ])
        cmd = [
            'asl_reg',
            f'-i {biascorr_name}',
            f'-s {t1_name}',
            f'--sbet {t1_brain_name}',
            f'-o {mask_dir}'
        ]
        subprocess.run(" ".join(cmd), shell=True)
        # apply transformation to the pve map
        for tissue in rois:
            roi_dir = mask_dir / roi
            create_dirs([roi_dir, ])
            pve_struct_name = fsl_anat_dir / PVE_NAMES[tissue]
            pve_asl_name = roi_dir / 'struct2asl.mat'
            struct2asl_name = roi_dir / 'struct2asl.mat'
            applyxfm(
                str(pve_struct_name),
                str(biascorr_name),
                str(struct2asl_name),
                str(pve_asl_name)
            )
            # threshold and binarise the ASL-space pve map
            mask_name = roi_dir / f'{tissue}_mask.nii.gz'
            fslmaths(str(pve_asl_name)).thr(0.9).bin().run(str(mask_name))
            # apply the mask to the calibration image
            masked_name = mask_dir / f'{tissue}_masked.nii.gz'
            fslmaths(str(biascorr_name)).mul(str(mask_name)).run(str(masked_name))
            # update dictionary
            new_dict = {
                f'{calib_name_stem}_{tissue}_masked': str(masked_name)
            }
            important_dict.update(new_dict)
    # save json
    with open(json_name, 'w') as fp:
        json.dump(important_dict, fp, sort_keys=True, indent=4)
    return (subject_dir, 0)