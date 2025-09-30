"""
Perform the setup necessary for estimating empiricial banding. This
includes finding necessary files, creating results directories
and running fsl_anat on the structural image.
"""

from functools import partial
import subprocess as sp
from pathlib import Path
import os.path as op

import nibabel as nb
import regtricks as rt
from fsl.wrappers import bet, fslmaths, fslroi

from hcpasl import distortion_correction
from hcpasl.bias_estimation import bias_estimation, register_fmap
from hcpasl.tissue_masks import generate_tissue_mask, generate_tissue_mask_in_ref_space
from hcpasl.utils import linear_asl_reg, load_asl_params, compute_split_indices


def setup(subject_dir):
    """
    Perform the initial set up for the empirical banding estimation pipeline.

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
    subject_dir = Path(subject_dir)
    subject_dir.resolve(strict=True)
    subid = subject_dir.parts[-1]

    # create directories
    calib_dirs = [subject_dir / f"ASL/calibration/{c}" for c in ("Calib0", "Calib1")]
    for d in calib_dirs:
        d.mkdir(exist_ok=True, parents=True)

    # mbPCASLhr_unproc and Structural_preproc directories
    mbpcasl_dir = subject_dir / "resources/mbPCASLhr_unproc/files"
    struct_dir = subject_dir / "resources/Structural_preproc/files"

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
    calib_dirs = [subject_dir / f"ASL/calibration/{c}" for c in ("Calib0", "Calib1")]
    for d in calib_dirs:
        d.mkdir(exist_ok=True, parents=True)
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

    # split mbpcasl sequence into its calibration images (always 2 calibration images)
    calib_names = [d / f"calib{n}.nii.gz" for n, d in enumerate(calib_dirs)]
    asl_params = load_asl_params(mbpcasl)
    img = nb.load(str(mbpcasl))
    total_vols = img.shape[3] if img.ndim == 4 else 1
    _, _, calib_idxs = compute_split_indices(
        total_vols,
        tail_discard_vols=asl_params.tail_discard_vols,
    )
    fslroi(str(mbpcasl), str(calib_names[0]), calib_idxs[0], 1)
    fslroi(str(mbpcasl), str(calib_names[1]), calib_idxs[1], 1)

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


def generate_sdc_warp(
    asl_vol0_brain,
    struct,
    struct_brain,
    asl_mask,
    wmmask,
    asl2struct,
    fmap,
    fmapmag,
    fmapmagbrain,
    distcorr_dir,
    interpolation=3,
):
    """
    Generate susceptibility distortion correction warp via asl_reg.

    Args:
        asl_vol0_brain: path to first volume of ASL series, brain-extracted
        struct: path to T1 image, ac_dc_restore
        struct_brain: path to brain-extracted T1 image, ac_dc_restore_brain
        asl_mask: path to brain mask in ASL space
        wmmask: path to WM mask in T1 space
        asl2struct: regtricks.Registration for asl to structural
        fmap: path to topup's field map in rad/s
        fmapmag: path to topup's field map magnitude
        fmapmagbrain: path to topup's field map magnitude, brain only
        distcorr_dir: path to directory in which to place output

    Returns:
        n/a, file 'asl2struct_warp.nii.gz' is created in output directory
    """

    a2s_fsl = op.join(distcorr_dir, "asl2struct.mat")
    asl2struct.save_fsl(a2s_fsl, asl_vol0_brain, struct)

    # get linear registration from fieldmaps to structural
    fmap_struct_dir = op.join(distcorr_dir, "fmap_struct_reg")
    bbr_fmap2struct = register_fmap(
        fmapmag, fmapmagbrain, struct, struct_brain, fmap_struct_dir, wmmask
    )

    # apply linear registration to fmap, fmapmag and fmapmagbrain
    bbr_fmap2struct = rt.Registration.from_flirt(
        bbr_fmap2struct, src=fmapmag, ref=struct
    )
    fmap_struct, fmapmag_struct, fmapmagbrain_struct = [
        op.join(fmap_struct_dir, f"fmap{ext}_struct.nii.gz")
        for ext in ("", "mag", "magbrain")
    ]
    for fmap_name, fmapstruct_name in zip(
        (fmap, fmapmag, fmapmagbrain),
        (fmap_struct, fmapmag_struct, fmapmagbrain_struct),
    ):
        fmapstruct_img = bbr_fmap2struct.apply_to_image(
            fmap_name, struct, order=interpolation
        )
        nb.save(fmapstruct_img, fmapstruct_name)

    # run asl_reg using pre-registered fieldmap images
    cmd = (
        "asl_reg -i {} -o {} ".format(asl_vol0_brain, distcorr_dir)
        + "-s {} --sbet={} -m {} ".format(struct, struct_brain, asl_mask)
        + "--tissseg={} --imat={} --finalonly ".format(wmmask, a2s_fsl)
        + "--fmap={} --fmapmag={} ".format(fmap_struct, fmapmag_struct)
        + "--fmapmagbrain={} --nofmapreg ".format(fmapmagbrain_struct)
        + "--echospacing=0.00057 --pedir=y"
    )
    sp.run(cmd)


def setup_empirical_estimation(
    subject_dir,
    coeffs_path,
    rois=["wm"],
    interpolation=3,
    ignore_dropouts=False,
    force_refresh=True,
):
    # find files and separate calibration images from mbPCASL sequence
    names_dict = setup(subject_dir)
    calib0_name, calib1_name = [names_dict[f"calib{n}_name"] for n in (0, 1)]
    calib_stems = [c.stem.split(".")[0] for c in (calib0_name, calib1_name)]

    # create results directory
    suf = "_ignoredropouts" if ignore_dropouts else ""
    results_dirs = [
        names_dict[f"calib{n}_dir"] / f"SEbased_t1mask{suf}" for n in (0, 1)
    ]
    for d in results_dirs:
        d.mkdir(exist_ok=True, parents=True)

    t1_name, t1_brain_name = [names_dict[name] for name in ("t1_name", "t1_brain_name")]

    # generate white matter mask in T1 space for use in BBRegistration
    wm_mask = names_dict["aslt1_dir"] / "wm_mask.nii.gz"
    if not wm_mask.exists() or force_refresh:
        nb.save(generate_tissue_mask(names_dict["aparc_aseg"], "wm"), wm_mask)

    # setup distortion correction results directories
    distcorr_dir = names_dict["calib0_dir"].parent / "distortion_correction"
    gdc_dir = distcorr_dir / "gradient_unwarp"
    topup_dir = distcorr_dir / "topup"
    calib_distcorr_dirs = [r / "distortion_correction" for r in results_dirs]
    for d in [distcorr_dir, gdc_dir, topup_dir, *calib_distcorr_dirs]:
        d.mkdir(exist_ok=True, parents=True)

    # gradient distortion correction
    gdc_warp = gdc_dir / "fullWarp_abs.nii.gz"
    if not gdc_warp.exists() or force_refresh:
        distortion_correction.generate_gdc_warp(
            names_dict["calib0_name"], coeffs_path, gdc_dir, interpolation
        )

    # topup
    fmap, fmapmag, fmapmagbrain = [
        topup_dir / f"fmap{ext}.nii.gz" for ext in ("", "mag", "magbrain")
    ]
    if not all([f.exists() for f in (fmap, fmapmag, fmapmagbrain)]) or force_refresh:
        pa_ap_sefms = topup_dir / "merged_sefms.nii.gz"
        distortion_correction.stack_fmaps(
            names_dict["pa_sefm"], names_dict["ap_sefm"], pa_ap_sefms
        )
        topup_params = topup_dir / "topup_params.txt"
        distortion_correction.generate_topup_params(topup_params)
        topup_config = "b02b0.cnf"
        distortion_correction.generate_fmaps(
            pa_ap_sefms, topup_params, topup_config, topup_dir, gdc_warp
        )

    # load gdc warp
    gdc_warp_reg = rt.NonLinearRegistration.from_fnirt(
        coefficients=gdc_warp,
        src=calib0_name,
        ref=calib0_name,
        intensity_correct=True,
    )
    # apply gdc and sdc to both calibration images
    for calib_name, results_dir in zip((calib0_name, calib1_name), calib_distcorr_dirs):
        # apply gdc to the calibration image
        calib_name_stem = calib_name.stem.split(".")[0]
        gdc_calib_name = results_dir / f"{calib_name_stem}_gdc.nii.gz"
        if not gdc_calib_name.exists() or force_refresh:
            gdc_calib = gdc_warp_reg.apply_to_image(
                calib_name, calib_name, order=interpolation
            )
            nb.save(gdc_calib, gdc_calib_name)

        # estimate initial registration via asl_reg
        asl_lin_reg = results_dir / "label_control_reg_linear"
        asl_lin_reg.mkdir(exist_ok=True, parents=True)
        asl2struct_lin = asl_lin_reg / "asl2struct.mat"
        if not asl2struct_lin.exists() or force_refresh:
            linear_asl_reg(gdc_calib_name, asl_lin_reg, t1_name, t1_brain_name, wm_mask)
        init_linear = rt.Registration.from_flirt(
            asl2struct_lin, src=calib_name, ref=t1_name
        )

        # run bet - should I instead get mask from T1 here?
        bet_results = results_dir / "bet"
        bet_mask = results_dir / "bet_mask.nii.gz"
        if not bet_mask.exists() or force_refresh:
            betted_calibration = bet(
                str(gdc_calib_name), str(bet_results), g=0.2, f=0.2, m=True
            )

        # get susceptibility distortion correction warps
        asl_nonlin_reg = results_dir / "label_control_reg_nonlinear"
        asl_nonlin_reg.mkdir(exist_ok=True, parents=True)
        struct2asl, asl2struct_warp = [
            asl_nonlin_reg / n for n in ("struct2asl.mat", "asl2struct_warp.nii.gz")
        ]
        if (
            not all([f.exists() for f in (asl2struct_warp, struct2asl)])
            or force_refresh
        ):
            distortion_correction.generate_sdc_warp(
                str(gdc_calib_name),
                str(t1_name),
                t1_brain_name,
                bet_mask,
                wm_mask,
                init_linear,
                fmap,
                fmapmag,
                fmapmagbrain,
                str(asl_nonlin_reg),
            )

        # chain gradient and susceptibility distortion correction warps together
        asl2struct_warp_reg = rt.NonLinearRegistration.from_fnirt(
            coefficients=asl2struct_warp,
            src=calib_name,
            ref=t1_name,
            intensity_correct=True,
        )
        struct2asl_reg = rt.Registration.from_flirt(
            struct2asl, src=t1_name, ref=calib_name
        )
        dc_warp = rt.chain(gdc_warp_reg, asl2struct_warp_reg, struct2asl_reg)
        dc_calib_name = results_dir / f"{calib_name_stem}_dc.nii.gz"
        if not dc_calib_name.exists() or force_refresh:
            dc_calib = dc_warp.apply_to_image(
                calib_name, calib_name, order=interpolation
            )
            nb.save(dc_calib, dc_calib_name)

        # estimate the bias field
        bias_name = results_dir / f"{calib_name_stem}_biasfield.nii.gz"
        sebased_dir = results_dir / "sebased"
        if not bias_name.exists() or force_refresh:
            sebased_dir.mkdir(exist_ok=True, parents=True)
            bias_field = bias_estimation(
                dc_calib_name,
                "sebased",
                results_dir=sebased_dir,
                t1_name=t1_name,
                t1_brain_name=t1_brain_name,
                aparc_aseg=names_dict["aparc_aseg"],
                fmapmag=fmapmag,
                fmapmagbrain=fmapmagbrain,
                interpolation=interpolation,
                force_refresh=force_refresh,
                wmseg_name=wm_mask,
                struct2asl=struct2asl,
            )
            nb.save(bias_field, bias_name)
        else:
            bias_field = nb.load(bias_name)

        # bias correct the distortion-corrected calibration image
        bc_calib_name = results_dir / f"{calib_name_stem}_bc.nii.gz"
        if not bc_calib_name.exists() or force_refresh:
            fslmaths(str(dc_calib_name)).div(bias_field).run(str(bc_calib_name))

        # create mask directories
        roi_dirs = [results_dir / "masks" / roi for roi in rois]
        for d in roi_dirs:
            d.mkdir(exist_ok=True, parents=True)
        # load Dropouts
        dropouts_inv = nb.load(sebased_dir / "Dropouts_inv.nii.gz")
        for roi, roi_dir in zip(rois, roi_dirs):
            # get tissue masks in calibration image space
            base_mask_call = partial(
                generate_tissue_mask_in_ref_space,
                aparc_aseg=names_dict["aparc_aseg"],
                ref_img=bc_calib_name,
                struct2ref=struct2asl,
                order=0,
            )
            tissues = ("gm", "wm") if roi == "combined" else (roi,)
            names = [roi_dir / f"{t}_mask.nii.gz" for t in tissues]
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
                masks = [
                    nb.nifti1.Nifti1Image(
                        m.get_fdata() * dropouts_inv.get_fdata(),
                        affine=dropouts_inv.affine,
                    )
                    for m in masks
                ]
                # save
                [nb.save(m, n) for m, n in zip(masks, names)]
            # apply tissue masks to bias- and distortion- corrected images
            calib_masked_names = [
                roi_dir / f"{calib_name_stem}_{t}_masked.nii.gz" for t in tissues
            ]
            if not all(c.exists() for c in calib_masked_names) or force_refresh:
                [
                    fslmaths(str(bc_calib_name)).mul(mask).run(str(name))
                    for mask, name in zip(masks, calib_masked_names)
                ]
    return (subject_dir, 1)
