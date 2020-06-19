"""
Prepare T1-asl gridded space
Estimate PVs in this space
Prepare ventricle mask in this space for final calibration 
"""

import os.path as op 
import sys 
import argparse
import multiprocessing as mp 
import glob 

import scipy
import numpy as np 
import toblerone as tob 
import regtricks as rt 
import nibabel as nib 

from hcpasl.extract_fs_pvs import extract_fs_pvs

def generate_ventricle_mask(aparc_aseg, t1_asl):

    ref_spc = rt.ImageSpace(t1_asl)

    # get ventricles mask from aparc+aseg image
    aseg = nib.load(aparc_aseg).get_data()
    vent_mask = np.logical_or(
        aseg == 43, # left ventricle 
        aseg == 4   # right ventricle 
    ) 

    # erosion in t1 space for safety 
    vent_mask = scipy.ndimage.morphology.binary_erosion(vent_mask)

    # Resample to target space, re-threshold
    output = rt.Registration.identity().apply_to_array(
                vent_mask, aparc_aseg, ref_spc)
    output = (output > 0.8)
    return output 


def estimate_pvs(t1_dir, t1_asl):
    """
    Generate partial volume estimates from freesurfer segmentations of the cortex
    and subcortical structures.
    Args: 
        t1w_dir: path to subject's T1w directory, containing a T1w scan (eg 
            acdc_dc_restore), aparc+aseg FS volumetric segmentation and 
            fsaverage_32k surface directory
        asl: path to ASL image, used for setting resolution of output 
        fileroot: path basename for output, will add suffix GM/WM/CSF
    """    

    # Load the t1 image, aparc+aseg and surfaces from their expected 
    # names and locations within t1w_dir 
    surf_dict = {}
    t1 = op.join(t1_dir, 'T1w_acpc_dc.nii.gz')
    aparc_aseg = op.join(t1_dir, 'aparc+aseg.nii.gz')
    for k,n in zip(['LWS', 'LPS', 'RPS', 'RWS'], 
                   ['L.white.32k_fs_LR.surf.gii',
                    'L.pial.32k_fs_LR.surf.gii', 
                    'R.pial.32k_fs_LR.surf.gii', 
                    'R.white.32k_fs_LR.surf.gii']):
        paths = glob.glob(op.join(t1_dir, 'fsaverage_LR32k', '*' + n))
        assert len(paths) == 1, f'Found multiple surfaces for {k}'
        surf_dict[k] = paths[0]

    # Generate a single 4D volume of PV estimates, stacked GM/WM/CSF
    pvs_stacked = extract_fs_pvs(aparc_aseg, surf_dict, t1_asl, 
        superfactor=2, cores=mp.cpu_count())

    return pvs_stacked


def main(study_dir, sub_id):

    # for debug, re-use intermediate results
    force_refresh = True 

    sub_base = op.abspath(op.join(study_dir, sub_id))
    t1_dir = op.join(sub_base, "T1w")
    asl = op.join(sub_base, "ASL", "TIs", "tis.nii.gz")
    struct = op.join(t1_dir, "T1w_acpc_dc_restore.nii.gz")

    # Create ASL-gridded version of T1 image 
    t1_asl_grid = op.join(t1_dir, "ASL", "reg", 
                          "ASL_grid_T1w_acpc_dc_restore.nii.gz")
    if not op.exists(t1_asl_grid) or force_refresh:
        asl_spc = rt.ImageSpace(asl)
        t1_spc = rt.ImageSpace(struct)
        t1_asl_grid_spc = t1_spc.resize_voxels(asl_spc.vox_size / t1_spc.vox_size)
        nib.save(
            rt.Registration.identity().apply_to_image(struct, t1_asl_grid_spc), 
            t1_asl_grid)

    # Create a ventricle CSF mask in T1 ASL space 
    ventricle_mask = op.join(sub_base, "T1w", "ASL", "PVEs",
                             "vent_csf_mask.nii.gz")
    if not op.exists(ventricle_mask) or force_refresh: 
        aparc_aseg = op.join(t1_dir, "aparc+aseg.nii.gz")
        vmask = generate_ventricle_mask(aparc_aseg, t1_asl_grid)
        rt.ImageSpace.save_like(t1_asl_grid, vmask, ventricle_mask)

    # Estimate PVs in T1 ASL space 
    pv_gm = op.join(sub_base, "T1w", "ASL", "PVEs", "pve_GM.nii.gz")
    if not op.exists(pv_gm) or force_refresh:
        aparc_seg = op.join(t1_dir, "aparc+aseg.nii.gz")
        pvs_stacked = estimate_pvs(t1_dir, t1_asl_grid)

        # Save output with tissue suffix 
        fileroot = op.join(sub_base, "T1w", "ASL", "PVEs", "pve")
        for idx, suffix in enumerate(['GM', 'WM', 'CSF']):
            p = "{}_{}.nii.gz".format(fileroot, suffix)
            rt.ImageSpace.save_like(t1_asl_grid, pvs_stacked.dataobj[...,idx], p)


if __name__ == "__main__":

    # study_dir = 'HCP_asl_min_req'
    # sub_number = 'HCA6002236'
    # sys.argv[1:] = ('%s %s' % (study_dir, sub_number)).split()

    # argument handling
    parser = argparse.ArgumentParser(description="Create T1-ASL space, estimate PVs, generate CSF ventricle mask")
    parser.add_argument(
        "study_dir",
        help="Path of the base study directory."
    )
    parser.add_argument(
        "sub_number",
        help="Subject number."
    )

    args = parser.parse_args()
    main(args.study_dir, args.sub_number)