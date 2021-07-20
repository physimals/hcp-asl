"""
Perform the setup necessary for estimating the MT Effect. This 
includes finding necessary files, creating results directories 
and running fsl_anat on the structural image.
"""

from pathlib import Path
from itertools import product
from functools import partial

from fsl.wrappers.fsl_anat import fsl_anat
from fsl.wrappers import fslmaths, LOAD, bet
from fsl.data.image import Image

import regtricks as rt
import nibabel as nb

from hcpasl import distortion_correction
from hcpasl.bias_estimation import bias_estimation, METHODS
from hcpasl.utils import (create_dirs, linear_asl_reg, setup,
                         binarise, get_ventricular_csf_mask)
from hcpasl.tissue_masks import (generate_tissue_mask, 
                                 generate_tissue_mask_in_ref_space)

def setup_mtestimation(subject_dir, coeffs_path, rois=['wm',],  
                        interpolation=3, ignore_dropouts=False,
                        force_refresh=True):

    try:
        # find files and separate calibration images from mbPCASL sequence
        names_dict = setup(subject_dir)
        calib0_name, calib1_name = [names_dict[f"calib{n}_name"] for n in (0, 1)]
        calib_stems = [c.stem.split(".")[0] for c in (calib0_name, calib1_name)]

        # create results directory
        suf = "_ignoredropouts" if ignore_dropouts else ""
        results_dirs = [names_dict[f'calib{n}_dir']/f"SEbased_MT_t1mask{suf}" for n in (0, 1)]
        create_dirs(results_dirs)

        t1_name, t1_brain_name = [names_dict[name] for name in ("t1_name", "t1_brain_name")]

        # generate white matter mask in T1 space for use in BBRegistration
        wm_mask = names_dict["aslt1_dir"]/"wm_mask.nii.gz"
        if not wm_mask.exists() or force_refresh:
            nb.save(generate_tissue_mask(names_dict["aparc_aseg"], "wm"), wm_mask)

        # setup distortion correction results directories
        distcorr_dir = names_dict["calib0_dir"].parent/"DistCorr"
        gdc_dir = distcorr_dir/"gradient_unwarp"
        topup_dir = distcorr_dir/"topup"
        calib_distcorr_dirs = [r/"DistCorr" for r in results_dirs]
        create_dirs((distcorr_dir, gdc_dir, topup_dir, *calib_distcorr_dirs))

        # gradient distortion correction
        gdc_warp = gdc_dir/"fullWarp_abs.nii.gz"
        if not gdc_warp.exists() or force_refresh:
            distortion_correction.generate_gdc_warp(
                names_dict["calib0_name"], coeffs_path, gdc_dir, interpolation
            )
        
        # topup
        fmap, fmapmag, fmapmagbrain = [
            topup_dir/f"fmap{ext}.nii.gz" for ext in ("", "mag", "magbrain")
        ]
        if not all([f.exists() for f in (fmap, fmapmag, fmapmagbrain)]) or force_refresh:
            pa_ap_sefms = topup_dir/"merged_sefms.nii.gz"
            distortion_correction.stack_fmaps(
                names_dict["pa_sefm"], names_dict["ap_sefm"], pa_ap_sefms
            )
            topup_params = topup_dir/"topup_params.txt"
            distortion_correction.generate_topup_params(topup_params)
            topup_config = "b02b0.cnf"
            distortion_correction.generate_fmaps(
                pa_ap_sefms, topup_params, topup_config, topup_dir, gdc_warp
            )
        
        # load gdc warp
        gdc_warp_reg = rt.NonLinearRegistration.from_fnirt(coefficients=str(gdc_warp),
                                                           src=str(calib0_name),
                                                           ref=str(calib0_name),
                                                           intensity_correct=True)
        # apply gdc and epidc to both calibration images
        for calib_name, results_dir in zip((calib0_name, calib1_name), calib_distcorr_dirs):
            # apply gdc to the calibration image
            calib_name_stem = calib_name.stem.split(".")[0]
            gdc_calib_name = results_dir/f"{calib_name_stem}_gdc.nii.gz"
            if not gdc_calib_name.exists() or force_refresh:
                gdc_calib = gdc_warp_reg.apply_to_image(
                    str(calib_name), str(calib_name), order=interpolation
                )
                nb.save(gdc_calib, gdc_calib_name)
            
            # estimate initial registration via asl_reg
            asl_lin_reg = results_dir/"asl_reg_linear"
            asl_lin_reg.mkdir(exist_ok=True)
            asl2struct_lin = asl_lin_reg/"asl2struct.mat"
            if not asl2struct_lin.exists() or force_refresh:
                linear_asl_reg(gdc_calib_name, asl_lin_reg, t1_name, t1_brain_name, wm_mask)
            init_linear = rt.Registration.from_flirt(str(asl2struct_lin),
                                                    src=str(calib_name), ref=str(t1_name))
            
            # run bet - should I instead get mask from T1 here?
            bet_results = results_dir/"bet"
            bet_mask = results_dir/"bet_mask.nii.gz"
            if not bet_mask.exists() or force_refresh:
                betted_m0 = bet(str(gdc_calib_name), str(bet_results), g=0.2, f=0.2, m=True)
            
            # get epi distortion correction warps
            asl_nonlin_reg = results_dir/"asl_reg_nonlinear"
            asl_nonlin_reg.mkdir(exist_ok=True)
            struct2asl, asl2struct_warp = [
                asl_nonlin_reg/n for n in ("struct2asl.mat", "asl2struct_warp.nii.gz")
            ]
            if not all([f.exists() for f in (asl2struct_warp, struct2asl)]) or force_refresh:
                distortion_correction.generate_epidc_warp(
                    str(gdc_calib_name), str(t1_name), t1_brain_name, bet_mask, wm_mask,
                    init_linear, fmap, fmapmag, fmapmagbrain, str(asl_nonlin_reg)
                )
            
            # chain gradient and epi distortion correction warps together
            asl2struct_warp_reg = rt.NonLinearRegistration.from_fnirt(coefficients=str(asl2struct_warp),
                                                                      src=str(calib_name),
                                                                      ref=str(t1_name),
                                                                      intensity_correct=True)
            struct2asl_reg = rt.Registration.from_flirt(
                str(struct2asl), src=str(t1_name), ref=str(calib_name)
            )
            dc_warp = rt.chain(gdc_warp_reg, asl2struct_warp_reg, struct2asl_reg)
            dc_calib_name = results_dir/f"{calib_name_stem}_dc.nii.gz"
            if not dc_calib_name.exists() or force_refresh:
                dc_calib = dc_warp.apply_to_image(
                    str(calib_name), str(calib_name), order=interpolation
                )
                nb.save(dc_calib, dc_calib_name)
            
            # estimate the bias field
            bias_name = results_dir/f"{calib_name_stem}_bias.nii.gz"
            sebased_dir = results_dir/"sebased"
            if not bias_name.exists() or force_refresh:
                sebased_dir.mkdir(exist_ok=True)
                bias_field = bias_estimation(
                    dc_calib_name, "sebased", results_dir=sebased_dir, t1_name=t1_name,
                    t1_brain_name=t1_brain_name, aparc_aseg=names_dict["aparc_aseg"], 
                    fmapmag=fmapmag, fmapmagbrain=fmapmagbrain, interpolation=interpolation,
                    force_refresh=force_refresh, wmseg_name=wm_mask, struct2asl=struct2asl
                )
                nb.save(bias_field, bias_name)
            else:
                bias_field = nb.load(bias_name)
            
            # bias correct the distortion-corrected calibration image
            bc_calib_name = results_dir/f"{calib_name_stem}_restore.nii.gz"
            if not bc_calib_name.exists() or force_refresh:
                fslmaths(str(dc_calib_name)).div(bias_field).run(str(bc_calib_name))

            # create mask directories
            roi_dirs = [results_dir/"masks"/roi for roi in rois]
            create_dirs(roi_dirs)
            # load Dropouts
            dropouts_inv = nb.load(sebased_dir/"Dropouts_inv.nii.gz")
            for roi, roi_dir in zip(rois, roi_dirs):
                # get tissue masks in calibration image space
                base_mask_call = partial(
                    generate_tissue_mask_in_ref_space,
                    aparc_aseg=names_dict["aparc_aseg"],
                    ref_img=bc_calib_name,
                    struct2ref=struct2asl,
                    order=0
                )
                tissues = ("gm", "wm") if roi=="combined" else (roi,)
                names = [roi_dir/f"{t}_mask.nii.gz" for t in tissues]
                if not all(n.exists() for n in names) or force_refresh:
                    if roi == "csf":
                        masks = [base_mask_call(tissue=t, erode=True) for t in tissues]
                    else:
                        masks = [base_mask_call(tissue=t) for t in tissues]
                    [nb.save(m, n) for m, n in zip(masks, names)]
                else:
                    masks = [nb.load(n) for n in names]
                if ignore_dropouts:
                    # ignore dropout voxels
                    masks = [nb.nifti1.Nifti1Image(
                                            m.get_fdata()*dropouts_inv.get_fdata(), 
                                            affine=dropouts_inv.affine)
                                for m in masks]
                    # save
                    [nb.save(m, n) for m, n in zip(masks, names)]
                # apply tissue masks to bias- and distortion- corrected images
                calib_masked_names = [roi_dir/f"{calib_name_stem}_{t}_masked.nii.gz"
                                    for t in tissues]
                if not all(c.exists() for c in calib_masked_names) or force_refresh:
                    [fslmaths(str(bc_calib_name)).mul(mask).run(str(name))
                                for mask, name in zip(masks, calib_masked_names)]
        return (subject_dir, 1)
    except Exception as e:
        print(e)
        return (subject_dir, e)
    