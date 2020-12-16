"""
This contains a range of functions required to correct for the 
Magnetisation Transfer effect visible in the HCP data.
"""

import json
from pathlib import Path
from fsl.wrappers import fslmaths, LOAD, bet, fast
from fsl.data.image import Image
import numpy as np
from .initial_bookkeeping import create_dirs
from .tissue_masks import generate_tissue_mask
from .distortion_correction import register_fmap, generate_asl_mask
import subprocess
import regtricks as rt
import nibabel as nb
import pandas as pd

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
    json_name = subject_dir / 'ASL/ASL.json'
    if json_name.exists():
        with open(json_name, 'r') as infile:
            json_dict = json.load(infile)
    else:
        raise Exception(f'File {json_name} does not exist.' + 
                        'Please run initial_processing() first.')
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
    with open(Path(old_dict['json_name']), 'w') as fp:
        json.dump(old_dict, fp, sort_keys=True, indent=4)

def parse_LUT(LUT_name):
    """
    Parse a FreeSurfer Lookup-Table returning the desired label 
    names.

    Parameters
    ----------
    LUT_name : str
        Path to the Lookup-Table
    
    Returns
    -------
    labels : list of ints
        List of int labels from the LUT
    """
    lut = pd.read_csv(LUT_name, header=None)
    labels = [int(row[0].split(" ")[0]) for _, row in lut[1::2].iterrows()]
    return labels

def correct_M0(subject_dir, mt_factors, wmparc, ribbon, 
               corticallut, subcorticallut, interpolation=3):
    """
    Correct the M0 images.
    
    For each of the subject's two calibration images:
    #. Use BET on the image;
    #. Use FAST on the brain-extracted image to obtain the bias-field;
    #. Perform bias correction;
    #. Multiply by the provided 'mt_factors' for MT-effect correction.
    
    Parameters
    ----------
    subject_dir : pathlib.Path
        Path to the subject's base directory.
    mt_factors : pathlib.Path
        Path to the empirically estimated MT correction 
        scaling factors.
    """
    # load json containing info on where files are stored
    json_dict = load_json(subject_dir/"hcp_asl")
    
    # do for both m0 images for the subject, calib0 and calib1
    calib_names = [json_dict['calib0_img'], json_dict['calib1_img']]

    # get structural image names
    struct_name, struct_brain_name = [json_dict[key] for key in ("T1w_acpc", "T1w_acpc_brain")]

    # generate white matter mask in T1w space for use in registration
    t1reg_dir = Path(json_dict["structasl"])/"reg"
    t1reg_dir.mkdir(exist_ok=True, parents=True)
    aparc_aseg = Path(json_dict["T1w_dir"])/"aparc+aseg.nii.gz"
    wmmask_img = generate_tissue_mask(aparc_aseg, "wm")
    wmmask_name = t1reg_dir/"wmmask.nii.gz"
    nb.save(wmmask_img, wmmask_name)

    # find gradient distortion correction warp and fieldmaps
    gdc_name = Path(json_dict['ASL_dir'])/'gradient_unwarp/fullWarp_abs.nii.gz'
    gdc_warp = rt.NonLinearRegistration.from_fnirt(
        str(gdc_name), calib_names[0], calib_names[0],
        intensity_correct=True, constrain_jac=(0.01, 100)
    )
    topup_dir = Path(json_dict["ASL_dir"])/"topup"
    fmap, fmapmag, fmapmagbrain = [topup_dir/f"fmap{ext}.nii.gz" 
                                   for ext in ('', 'mag', 'magbrain')]

    # register fieldmaps to structural image for use with epi_reg later
    fmap_struct_dir = topup_dir/"fmap_struct_reg"
    bbr_fmap2struct = register_fmap(fmapmag, fmapmagbrain, struct_name, 
                                    struct_brain_name, fmap_struct_dir, 
                                    wmmask_name)
    bbr_fmap2struct = rt.Registration.from_flirt(bbr_fmap2struct, 
                                                 src=str(fmapmag), 
                                                 ref=str(struct_name))
    fmap_struct, fmapmag_struct, fmapmagbrain_struct = [
        fmap_struct_dir/f"fmap{ext}_struct.nii.gz" for ext in ("", "mag", "magbrain")
    ]
    for fmap_name, fmapstruct_name in zip((fmap, fmapmag, fmapmagbrain),
                                          (fmap_struct, fmapmag_struct, fmapmagbrain_struct)):
        fmapstruct_img = bbr_fmap2struct.apply_to_image(str(fmap_name), 
                                                        str(struct_name), 
                                                        order=interpolation)
        nb.save(fmapstruct_img, fmapstruct_name)    

    for calib_name in calib_names:
        # get calib_dir and other info
        calib_path = Path(calib_name)
        calib_dir = calib_path.parent
        calib_name_stem = calib_path.stem.split('.')[0]

        # apply gradient distortion correction to the calibration image
        gdc_calib_img = gdc_warp.apply_to_image(calib_name, calib_name, order=interpolation)
        distcorr_dir = calib_dir/"DistCorr"
        distcorr_dir.mkdir(exist_ok=True)
        gdc_calib_name = distcorr_dir/f"gdc_{calib_name_stem}.nii.gz"
        nb.save(gdc_calib_img, gdc_calib_name)

        # apply mt scaling factors to the gradient distortion-corrected calibration image
        mt_sfs = np.loadtxt(mt_factors)
        assert (len(mt_sfs) == gdc_calib_img.shape[2])
        mt_gdc_calib_img = nb.nifti1.Nifti1Image(gdc_calib_img.get_fdata()*mt_sfs, 
                                                 gdc_calib_img.affine)
        mtcorr_dir = calib_dir/"MTCorr"
        mtcorr_dir.mkdir(exist_ok=True)
        mt_gdc_calib_name = mtcorr_dir/f"mtcorr_gdc_{calib_name_stem}.nii.gz"
        nb.save(mt_gdc_calib_img, mt_gdc_calib_name)

        # get registration to structural
        # initial
        asl_reg_cmd = [
            "asl_reg",
            "-i", mt_gdc_calib_name, "-o", distcorr_dir,
            "-s", struct_name, f"--sbet={struct_brain_name}",
            f"--tissseg={wmmask_name}", "--mainonly"
        ]
        print("Running first asl_reg")
        subprocess.run(asl_reg_cmd, check=True)

        # get brain mask in calibration image space
        calib2struct_init = distcorr_dir/"asl2struct.mat"
        struct2calib_reg = rt.Registration.from_flirt(str(calib2struct_init),
                                                      str(calib_name),
                                                      str(struct_name)
                                                      ).inverse()
        mask = generate_asl_mask(str(struct_brain_name), 
                                 str(calib_name), 
                                 struct2calib_reg.inverse())
        mask_name = calib_dir/"mask.nii.gz"
        rt.ImageSpace.save_like(calib_name, mask, str(mask_name))

        # refined
        asl_reg_cmd = [
            "asl_reg",
            "-i", str(mt_gdc_calib_name), "-o", str(distcorr_dir),
            "-s", str(struct_name), f"--sbet={str(struct_brain_name)}",
            f"--tissseg={str(wmmask_name)}", "--echospacing=0.00057", "--pedir=y",
            f"--fmap={str(fmap_struct)}", f"--fmapmag={str(fmapmag_struct)}", 
            f"--fmapmagbrain={str(fmapmagbrain_struct)}", "--nofmapreg", 
            f"--imat={str(calib2struct_init)}", "--finalonly", "-m", str(mask_name)
        ]
        print("Running second asl_reg")
        subprocess.run(asl_reg_cmd, check=True)
        struct2calib_name = distcorr_dir/"struct2asl.mat"
        struct2calib_reg = rt.Registration.from_flirt(str(struct2calib_name), 
                                                      str(struct_name),
                                                      str(calib_name))
        calib2struct_name = distcorr_dir/"asl2struct_warp.nii.gz"
        calib2struct_warp = rt.NonLinearRegistration.from_fnirt(str(calib2struct_name),
                                                                str(calib_name),
                                                                str(struct_name),
                                                                intensity_correct=True,
                                                                constrain_jac=(0.01, 100))
        
        # register fmapmag to calibration image space
        fmap2calib_reg = rt.chain(bbr_fmap2struct, struct2calib_reg)
        fmapmag_calibspc = fmap2calib_reg.apply_to_image(str(fmapmag),
                                                         str(calib_name),
                                                         order=interpolation)
        biascorr_dir = calib_dir/"BiasCorr"
        sebased_dir = biascorr_dir/"SEbased"
        sebased_dir.mkdir(parents=True, exist_ok=True)
        fmapmag_cspc_name = sebased_dir/f"fmapmag_{calib_name_stem}spc.nii.gz"
        nb.save(fmapmag_calibspc, fmapmag_cspc_name)

        # apply all distortion corrections to calibration image
        dc_warp = rt.chain(gdc_warp, calib2struct_warp, struct2calib_reg)
        dc_calib = dc_warp.apply_to_image(calib_name, 
                                          calib_name,
                                          order=interpolation)
        dc_calib_name = distcorr_dir/f"dc_{calib_name_stem}.nii.gz"
        nb.save(dc_calib, dc_calib_name)

        # get brain mask in calibration image space
        fs_brainmask = Path(json_dict["T1w_dir"])/"brainmask_fs.nii.gz"
        aslfs_mask_name = calib_dir/"aslfs_mask.nii.gz"
        aslfs_mask = struct2calib_reg.apply_to_image(src=str(fs_brainmask), 
                                                     ref=str(calib_name),
                                                     order=0)
        aslfs_mask = nb.nifti1.Nifti1Image(np.where(aslfs_mask.get_fdata()>0., 1., 0.),
                                           affine=gdc_calib_img.affine)
        nb.save(aslfs_mask, aslfs_mask_name)

        # get sebased bias estimate
        sebased_cmd = [
            "get_sebased_bias",
            "-i", dc_calib_name, "-f", fmapmag_cspc_name,
            "-m", aslfs_mask_name, "-o", sebased_dir, 
            "--ribbon", ribbon, "--wmparc", wmparc,
            "--corticallut", corticallut, "--subcorticallut", subcorticallut,
            "--struct2calib", struct2calib_name, "--structural", struct_name, 
            "--debug"
        ]
        subprocess.run(sebased_cmd, check=True)

        # apply dilall to bias estimate
        bias_name = sebased_dir/"sebased_bias_dil.nii.gz"
        dilall_name = biascorr_dir/f"{calib_name_stem}_bias.nii.gz"
        dilall_cmd = ["fslmaths", bias_name, "-dilall", dilall_name]
        subprocess.run(dilall_cmd, check=True)

        # apply mt scaling factors to dc calib bias field
        bias_img = nb.load(dilall_name)
        bc_calib = nb.nifti1.Nifti1Image(dc_calib.get_fdata() / bias_img.get_fdata(),
                                         dc_calib.affine)
        biascorr_name = biascorr_dir / f'{calib_name_stem}_restore.nii.gz'
        mt_bc_calib = nb.nifti1.Nifti1Image(bc_calib.get_fdata()*mt_sfs,
                                            bc_calib.affine)
        mtcorr_name = mtcorr_dir / f'{calib_name_stem}_mtcorr.nii.gz'
        [nb.save(image, name) for image, name in zip((bc_calib, mt_bc_calib),
                                                     (biascorr_name, mtcorr_name))]
        
        # add locations of above files to the json
        important_names = {
            f'{calib_name_stem}_bias' : str(dilall_name),
            f'{calib_name_stem}_bc' : str(biascorr_name),
            f'{calib_name_stem}_mc' : str(mtcorr_name)
        }
        update_json(important_names, json_dict)