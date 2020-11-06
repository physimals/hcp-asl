"""
Perform the setup necessary for estimating the MT Effect. This 
includes finding necessary files, creating results directories 
and running fsl_anat on the structural image.
"""

from pathlib import Path
import os
from fsl.wrappers.misc import fslroi
from fsl.wrappers.fsl_anat import fsl_anat
from fsl.wrappers.flirt import applyxfm
from fsl.wrappers.fnirt import applywarp, invwarp
from fsl.wrappers import fslmaths, LOAD, bet, fast
from fsl.data.image import Image
from fsl.data import atlases
import json
from hcpasl.initial_bookkeeping import create_dirs
import subprocess

PVE_NAMES = {
    'csf': 'T1_fast_pve_0.nii.gz',
    'gm': 'T1_fast_pve_1.nii.gz',
    'wm': 'T1_fast_pve_2.nii.gz'
}
PVE_THRESHOLDS = {
    'csf': 0.9,
    'gm': 0.7,
    'wm': 0.9,
    'combined': 0.7
}

def setup_mtestimation(subject_dir, rois=['wm',], biascorr_method='calib', force_refresh=False):
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
    if not fsl_anat_dir.with_suffix('.anat').exists():
        print(f"{fsl_anat_dir.with_suffix('.anat')} exists: {fsl_anat_dir.exists()}")
        fsl_anat(str(struc_name), str(fsl_anat_dir), clobber=True, nosubcortseg=True)
    fsl_anat_dir = fsl_anat_dir.parent / f'{fsl_anat_dir.stem}.anat'
    t1_name = fsl_anat_dir / 'T1_biascorr.nii.gz'
    t1_brain_name = fsl_anat_dir / 'T1_biascorr_brain.nii.gz'

    # get ventricles mask if csf
    if 'csf' in rois:
        vent_t1_name = fsl_anat_dir / 'ventricles_mask.nii.gz'
        if (not vent_t1_name.exists()) or force_refresh:
            # initialise atlas list
            atlases.rescanAtlases()
            harv_ox_prob_2mm = atlases.loadAtlas(
                'harvardoxford-subcortical', 
                resolution=2.0
            )
            vent_img = Image(
                harv_ox_prob_2mm.data[:, :, :, 2] + harv_ox_prob_2mm.data[:, :, :, 13],
                header=harv_ox_prob_2mm.header
            )
            vent_img = fslmaths(vent_img).thr(0.1).bin().ero().run(LOAD)
            # we already have registration from T1 to MNI
            struc2mni_coeff = fsl_anat_dir / 'T1_to_MNI_nonlin_coeff.nii.gz'
            mni2struc_coeff = fsl_anat_dir / 'MNI_to_T1_nonlin_coeff.nii.gz'
            invwarp(str(struc2mni_coeff), str(t1_brain_name), str(mni2struc_coeff))
            # apply warp to ventricles image
            applywarp(vent_img, str(t1_brain_name), str(vent_t1_name), warp=str(mni2struc_coeff))
            # re-threshold
            vent_t1 = fslmaths(str(vent_t1_name)).thr(PVE_THRESHOLDS['csf']).bin().run(LOAD)
            # mask pve estimate by ventricles mask
            fslmaths(str(fsl_anat_dir / PVE_NAMES['csf'])).mas(vent_t1).run(str(vent_t1_name))

    # bias-field correction
    for calib_name in (calib0_name, calib1_name):
        calib_name_stem = calib_name.stem.split('.')[0]
        # run bet
        betted_m0 = bet(str(calib_name), LOAD, g=0.2, f=0.2, m=True)
        # create directories to save results
        biascorr_dir = calib_name.parent / 'BiasCorr'
        if biascorr_method == 'calib':
            fast_dir = calib_name.parent / 'Fast'
            create_dirs([biascorr_dir, fast_dir])
        else:
            create_dirs([biascorr_dir,])
        # obtain registration from structural to calibration image
        mask_dir = calib_name.parent / 'masks'
        create_dirs([mask_dir, ])
        struct2asl_name = (mask_dir/f'struct2asl_{biascorr_method}.mat')
        if (not struct2asl_name.exists()) or force_refresh:
            struct2asl_name_temp = mask_dir/'struct2asl.mat'
            cmd = [
                'asl_reg',
                f'-i {calib_name}',
                f'-s {t1_name}',
                f'--sbet {t1_brain_name}',
                f'-o {mask_dir}'
            ]
            subprocess.run(" ".join(cmd), shell=True)
            os.replace(struct2asl_name_temp, struct2asl_name)
        # perform bias correction
        if biascorr_method == 'calib':
            biascorr_name = biascorr_dir / f'{calib_name_stem}_restore.nii.gz'
        else:
            biascorr_name = biascorr_dir / f'T1_biascorr_{calib_name_stem}.nii.gz'
        if (not biascorr_name.exists()) or force_refresh:
            if biascorr_method == 'calib':
                fast_base = fast_dir / 'fast'
                fast_results = fast(
                    betted_m0['output'], # output of bet
                    out=str(fast_base), 
                    type=3, # image type, 3=PD image
                    b=True, # output estimated bias field
                    nopve=True # don't need pv estimates
                )
                calib_bias_name = fast_base / 'fast_bias.nii.gz'
            else:
                # get bias field from T1 image's fsl_anat
                bias_name = fsl_anat_dir / 'T1_fast_bias.nii.gz'
                # apply registration to the T1 bias field
                calib_bias_name = biascorr_dir / 'T1_bias_field_calibspc.nii.gz'
                applyxfm(
                    str(bias_name),
                    str(calib_name),
                    str(struct2asl_name),
                    str(calib_bias_name)
                )
            fslmaths(str(calib_name)).div(str(calib_bias_name)).run(str(biascorr_name))
        # apply transformation to the pve map
        for tissue in rois:
            roi_dir = mask_dir / tissue
            create_dirs([roi_dir, ])
            if tissue == 'combined':
                # check that gm and wm pves already exist
                gm_pve = mask_dir / 'gm/pve_gm.nii.gz'
                wm_pve = mask_dir / 'wm/pve_wm.nii.gz'
                if not (gm_pve.exists() and wm_pve.exists()):
                    raise Exception("WM and GM PVEs don't exist.")
                gm_mask_name = roi_dir / f'gm_mask_{biascorr_method}.nii.gz'
                wm_mask_name = roi_dir / f'wm_mask_{biascorr_method}.nii.gz'
                fslmaths(str(gm_pve)).thr(PVE_THRESHOLDS[tissue]).bin().run(str(gm_mask_name))
                fslmaths(str(wm_pve)).thr(PVE_THRESHOLDS[tissue]).bin().run(str(wm_mask_name))
                gm_masked_name = roi_dir / f'{calib_name_stem}_gm_masked_{biascorr_method}'
                wm_masked_name = roi_dir / f'{calib_name_stem}_wm_masked_{biascorr_method}'
                fslmaths(str(biascorr_name)).mul(str(gm_mask_name)).mul(betted_m0['output_mask']).run(str(gm_masked_name))
                fslmaths(str(biascorr_name)).mul(str(wm_mask_name)).mul(betted_m0['output_mask']).run(str(wm_masked_name))
                # update dictionary
                new_dict = {
                    f'{calib_name_stem}_{tissue}_masked': [
                        str(gm_masked_name), 
                        str(wm_masked_name)
                    ]
                }
            else:
                if tissue == 'csf':
                    pve_struct_name = vent_t1_name
                else:
                    pve_struct_name = fsl_anat_dir / PVE_NAMES[tissue]
                pve_asl_name = roi_dir / f'pve_{tissue}.nii.gz'
                applyxfm(
                    str(pve_struct_name),
                    str(calib_name),
                    str(struct2asl_name),
                    str(pve_asl_name)
                )
                # threshold and binarise the ASL-space pve map
                mask_name = roi_dir / f'{tissue}_mask_{biascorr_method}.nii.gz'
                fslmaths(str(pve_asl_name)).thr(PVE_THRESHOLDS[tissue]).bin().run(str(mask_name))
                # apply the mask to the calibration image
                masked_name = mask_dir / f'{tissue}_masked_{biascorr_method}.nii.gz'
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