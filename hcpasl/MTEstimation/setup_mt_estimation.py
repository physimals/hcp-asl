"""
Perform the setup necessary for estimating the MT Effect. This 
includes finding necessary files, creating results directories 
and running fsl_anat on the structural image.
"""

from pathlib import Path
from itertools import product

from fsl.wrappers.fsl_anat import fsl_anat
from fsl.wrappers import fslmaths, LOAD, bet
from fsl.data.image import Image

import regtricks as rt
import nibabel as nb

from hcpasl.initial_bookkeeping import create_dirs
from hcpasl import distortion_correction
from hcpasl.bias_estimation import bias_estimation, METHODS
from hcpasl.utils import (create_dirs, linear_asl_reg, setup,
                         binarise, get_ventricular_csf_mask)

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

def setup_mtestimation(subject_dir, rois=['wm',], biascorr_method='calib', 
                        distcorr=True, coeffs_path=None, interpolation=3, 
                        force_refresh=True, fslanat_refresh=False):

    # make sure chose bias estimation method is supported
    assert (biascorr_method in METHODS), f"{biascorr_method} is not a supported method."

    # find files and separate calibration images from mbPCASL sequence
    names_dict = setup(subject_dir)
    calib0_name, calib1_name = [names_dict[f"calib{n}_name"] for n in (0, 1)]
    calib_stems = [c.stem.split(".")[0] for c in (calib0_name, calib1_name)]

    # create results directory
    suffix = "_newdistcorr" if distcorr else ""
    results_dirs = [names_dict[f'calib{n}_dir']/f"{biascorr_method}{suffix}" for n in (0, 1)]
    create_dirs(results_dirs)

    # run fsl_anat
    fslanatdir = names_dict['aslt1_dir']/"struc"
    if not fslanatdir.with_suffix('.anat').exists() or fslanat_refresh:
        fsl_anat(str(names_dict['t1_name']), str(fslanatdir), clobber=True, nosubcortseg=True)
    fslanatdir = fslanatdir.with_suffix(".anat")
    t1_name = fslanatdir/"T1_biascorr.nii.gz"
    t1_brain_name = fslanatdir/"T1_biascorr_brain.nii.gz"
    wmpve_name = fslanatdir/"T1_fast_pve_2.nii.gz"

    # perform distortion correction
    if distcorr:
        distcorr_dir = names_dict['calib0_dir'].parent/"DistCorr"
        calib_distcorr_dirs = [r/"DistCorr" for r in results_dirs]
        create_dirs((distcorr_dir, *calib_distcorr_dirs))
        # gradient distortion correction
        coeffs_path = Path(coeffs_path).resolve(strict=True)
        gdc_dir = distcorr_dir/"gradient_unwarp"
        gdc_dir.mkdir(exist_ok=True)
        gdc_warp = gdc_dir/"fullWarp_abs.nii.gz"
        if not gdc_warp.exists() or force_refresh:
            # run gradient_unwarp
            distortion_correction.generate_gdc_warp(
                names_dict['calib0_name'], coeffs_path, gdc_dir, interpolation
            )
        # epi distortion correction
        # stack images together for use with topup
        topup_dir = distcorr_dir/"topup"
        topup_dir.mkdir(exist_ok=True)
        pa_ap_sefms = topup_dir/"merged_sefms.nii.gz"
        if not pa_ap_sefms.exists() or force_refresh:
            distortion_correction.stack_fmaps(
                names_dict["pa_sefm"], names_dict["ap_sefm"], pa_ap_sefms
            )
        # generate topup params
        topup_params = topup_dir/"topup_params.txt"
        distortion_correction.generate_topup_params(topup_params)
        # run topup
        topup_config = "b02b0.cnf"
        fmap, fmapmag, fmapmagbrain = [
            topup_dir/f"fmap{ext}.nii.gz" for ext in ('', 'mag', 'magbrain')
        ]
        if not all([f.exists() for f in (fmap, fmapmag, fmapmagbrain)]) or force_refresh:
            distortion_correction.generate_fmaps(
                pa_ap_sefms, topup_params, topup_config, topup_dir
            )
        # apply gdc and epi dc to both calibration images
        gdc_warp_reg = rt.NonLinearRegistration.from_fnirt(
            str(gdc_warp), src=str(calib0_name), ref=str(calib0_name),
            intensity_correct=True, constrain_jac=(0.01, 100)
        )
        for calib_name, results_dir in zip((calib0_name, calib1_name), calib_distcorr_dirs):
            calib_name_stem = calib_name.stem.split('.')[0]
            # apply gdc to the calibration image
            gdc_calib_name = results_dir/f"{calib_name_stem}_gdc.nii.gz"
            if not gdc_calib_name.exists() or force_refresh:
                gdc_calib = gdc_warp_reg.apply_to_image(
                    str(calib_name), str(calib_name), order=interpolation
                )
                nb.save(gdc_calib, gdc_calib_name)
            # run asl_reg on the gdc calibration image
            asl_lin_reg = results_dir/"asl_reg_linear"
            asl_lin_reg.mkdir(exist_ok=True)
            asl2struct_lin = asl_lin_reg/"asl2struct.mat"
            if not asl2struct_lin.exists() or force_refresh:
                linear_asl_reg(gdc_calib_name, asl_lin_reg, t1_name, t1_brain_name, wmpve_name)
            init_linear = rt.Registration.from_flirt(str(asl2struct_lin),
                                                    src=str(calib_name), ref=str(t1_name))
            wmseg_name = asl_lin_reg/"wm_seg.nii.gz"
            # run bet
            bet_results = results_dir/"bet"
            bet_mask = results_dir/"bet_mask.nii.gz"
            if not bet_mask.exists() or force_refresh:
                betted_m0 = bet(str(gdc_calib_name), str(bet_results), g=0.2, f=0.2, m=True)
            # run together with asl_reg to get epi distortion correction warps
            asl_nonlin_reg = results_dir/"asl_reg_nonlinear"
            asl_nonlin_reg.mkdir(exist_ok=True)
            asl2struct_warp = asl_nonlin_reg/"asl2struct_warp.nii.gz"
            struct2asl = asl_nonlin_reg/"struct2asl.mat"
            if not all([f.exists() for f in (asl2struct_warp, struct2asl)]) or force_refresh:
                distortion_correction.generate_epidc_warp(
                    str(gdc_calib_name), str(t1_name), t1_brain_name, bet_mask, wmseg_name, init_linear, 
                    fmap, fmapmag, fmapmagbrain, str(asl_nonlin_reg)
                )
            # chain gradient and epi distortion correction warps together
            asl2struct_warp_reg = rt.NonLinearRegistration.from_fnirt(
                str(asl2struct_warp), src=str(calib_name), ref=str(t1_name), 
                intensity_correct=True, constrain_jac=(0.01, 100)
            )
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
        # set calib names to dc_calib so bias field is estimated 
        # on the distortion corrected images
        calib0_name, calib1_name = [
            distcorr_dir/f"calib{n}_dc.nii.gz" for n, distcorr_dir in enumerate(calib_distcorr_dirs)
        ]

    # estimate the bias field
    bias_names = [r/f"{biascorr_method}_bias.nii.gz" for r in results_dirs]
    if not all([b.exists() for b in bias_names]) or force_refresh:
        if biascorr_method == "calib":
            bias_fields = [
                bias_estimation(c, biascorr_method) for c in (calib0_name, calib1_name)
            ]
        elif biascorr_method == "t1":
            if distcorr:
                struct2asl_names = [
                    c/"asl_reg_nonlinear/struct2asl.mat" for c in calib_distcorr_dirs
                ]
            else:
                asl_reg_lin_dirs = [c/"asl_reg_lin" for c in results_dirs]
                create_dirs(asl_reg_lin_dirs)
                [linear_asl_reg(
                    calib_name=c, results_dir=a, t1_name=str(t1_name), 
                    t1_brain_name=str(t1_brain_name), wmpve_name=str(wmpve_name)
                ) for c, a in zip((calib0_name, calib1_name), asl_reg_lin_dirs)]
                struct2asl_names = [a/"struct2asl.mat" for a in asl_reg_lin_dirs]
            bias_fields = [
                bias_estimation(
                    c, biascorr_method, fslanatdir=fslanatdir, 
                    struct2asl=s2a, interpolation=interpolation
                ) 
                for c, s2a in zip((calib0_name, calib1_name), struct2asl_names)
            ]
        elif biascorr_method == "sebased":
            assert distcorr, "Should be performing distortion correction for sebased bias correction."
            struct2asl_names = [
                c/"asl_reg_nonlinear/struct2asl.mat" for c in calib_distcorr_dirs
            ]
            sebased_dirs = [r/"sebased" for r in results_dirs]
            create_dirs(sebased_dirs)
            bias_fields = [
                bias_estimation(
                    c, biascorr_method, results_dir=r, fslanatdir=fslanatdir, struct2asl=s,
                    fmapmag=fmapmag, fmapmagbrain=fmapmagbrain, wmseg_name=wmseg_name, 
                    interpolation=interpolation, force_refresh=force_refresh
                )
                for c, r, s in zip((calib0_name, calib1_name), sebased_dirs, struct2asl_names)
            ]
        # save bias field
        [nb.save(b, n) for b, n in zip(bias_fields, bias_names)]
    else:
        bias_fields = [nb.load(b) for b in bias_names]
    # use bias field to correct calibration images
    bc_calib_names = [r/f"{c}_restore.nii.gz" for r, c in zip(results_dirs, calib_stems)]
    if not all([bc_calib.exists() for bc_calib in bc_calib_names]) or force_refresh:
        bc_calib_imgs = [
            fslmaths(str(c)).div(b).run(str(bc_calib))
            for c, b, bc_calib in zip((calib0_name, calib1_name), bias_fields, bc_calib_names)
        ]

    # register T1 fast PVEs to the calibration image
    # if biascorr_method is "calib", struct2asl_names doesn't exist yet
    if biascorr_method == "calib":
        asl_reg_lin_dirs = [c/"asl_reg_lin" for c in results_dirs]
        create_dirs(asl_reg_lin_dirs)
        struct2asl_names = [a/"struct2asl.mat" for a in asl_reg_lin_dirs]
        if not all([s.exists() for s in struct2asl_names]) or force_refresh:
            [linear_asl_reg(
                calib_name=c, results_dir=a, t1_name=str(t1_name), 
                t1_brain_name=str(t1_brain_name), wmpve_name=str(wmpve_name)
            ) for c, a in zip((calib0_name, calib1_name), asl_reg_lin_dirs)]

    # load struct2asl registrations for both calibration images
    struct2asl_regs = [
        rt.Registration.from_flirt(str(s), str(t1_name), str(c)) 
        for s, c in zip(struct2asl_names, (calib0_name, calib1_name))
    ]
    # create roi results directories
    roi_dirs = [r/"masks"/roi for r, roi in product(results_dirs, rois)]
    create_dirs(roi_dirs)
    # load T1 space PVEs and apply registration to calibration images, 
    # masking with ventricles mask if ROI=CSF
    for roi in rois:
        tissues = ("gm", "wm") if roi=="combined" else (roi, )
        if roi == "csf":
            vent_t1_names = [r/"masks"/roi/"ventricles_mask_t1.nii.gz" for r in results_dirs]
            if not all([vent_name.exists() for vent_name in vent_t1_names]) or force_refresh:
                vent_t1_img = get_ventricular_csf_mask(fslanatdir, interpolation=interpolation)
                [vent_t1_img.save(str(vent_name)) for vent_name in vent_t1_names]

        calib_mask_names = [r/"masks"/roi/f"{t}_mask.nii.gz" for r, t in product(results_dirs, tissues)]
        if not all([c.exists() for c in calib_mask_names]) or force_refresh:
            # load tissue T1 space PVEs
            t1_pves = [Image(str(fslanatdir/PVE_NAMES[t])) for t in tissues]
            # mask by ventricles mask if roi == csf
            if roi == "csf":
                t1_pves = [
                    fslmaths(t1_pves[0]).mas(str(vent_t1)).run(LOAD) 
                    for vent_t1 in vent_t1_names
                ]
            # apply registrations to calibration space
            calib_pves = [
                s.apply_to_image(src=p, ref=str(c), order=interpolation)
                for (s, c), p in product(zip(struct2asl_regs, (calib0_name, calib1_name)), t1_pves)
            ]
            # save calibration space PVEs
            calib_pve_names = [r/"masks"/roi/f"{t}_pve.nii.gz" for r, t in product(results_dirs, tissues)]
            [c.save(str(n)) for c, n in zip(calib_pves, calib_pve_names)]
            # binarise PVEs to create calibration space tissue masks
            [
                binarise(p, threshold=PVE_THRESHOLDS[roi]).save(str(n)) 
                for n, p in zip(calib_mask_names, calib_pve_names)
            ]
            # apply mask to bias-corrected calibration image
            calib_masked_names = [
                r/"masks"/roi/f"{stem}_{t}_masked.nii.gz" 
                for (r, stem), t in product(zip(results_dirs, calib_stems), tissues)
            ]
            [
                fslmaths(str(calib_name)).mul(str(mask_name)).run(str(masked_name))
                for calib_name, (mask_name, masked_name) in product(
                    bc_calib_names, zip(calib_mask_names, calib_masked_names)
                )
            ]
