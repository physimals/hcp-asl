import subprocess
from pathlib import Path

import nibabel as nb
import numpy as np
import regtricks as rt
from fsl.wrappers import LOAD, bet, fast

from hcpasl.distortion_correction import register_fmap
from hcpasl.tissue_masks import generate_tissue_mask_in_ref_space


def bias_estimation_calib(calib_name):
    """
    Estimate the bias field from a calibration image using FAST.

    Parameters
    ----------
    calib_name: pathlib.Path
        Path to the calibration image from which to estimate the
        bias field.

    Returns
    -------
    bias_field: Nifti1Image
    """
    # perform brain extraction using BET
    betted_m0 = bet(str(calib_name), LOAD, g=0.2, f=0.2, m=True)
    # run FAST on BET-ed calibration image
    fast_results = fast(betted_m0["output"], out=LOAD, type=3, b=True, nopve=True)
    # extract and return bias field
    bias_field = fast_results["out_bias"]
    return bias_field


def bias_estimation_t1(calib_name, fslanatdir, struct2asl, interpolation=3):
    """
    Obtain a bias field estimate in calibration space by
    registering that obtained from the T1 image.

    Parameters
    ----------
    calib_name: pathlib.Path
        Path to the calibration image.
    fslanatdir: pathlib.Path
        Path to the T1 image's fslanat output.
    struct2asl: pathlib.Path
        Path to the .mat file which will be used to register
        the T1-estimated bias field to calibration image space.
    interpolation: int, default=3
        Order of interpolation to be used when registering
        the T1 image's bias field to the calibration image.

    Returns
    -------
    bias_field: Nifti1Image
    """
    # get t1 and bias field
    t1_name = fslanatdir / "T1_biascorr.nii.gz"
    t1_bias = fslanatdir / "T1_fast_bias.nii.gz"
    # load struct2asl registration
    struct2asl_reg = rt.Registration.from_flirt(struct2asl, src=t1_name, ref=calib_name)
    # apply registration to bias field to get it in calib space
    bias_field = struct2asl_reg.apply_to_image(t1_bias, calib_name, order=interpolation)
    return bias_field


def bias_estimation_sebased(
    calib_name,
    struct2asl,
    wmseg_name,
    results_dir,
    t1_name,
    t1_brain_name,
    aparc_aseg,
    fmapmag,
    fmapmagbrain,
    interpolation=3,
    force_refresh=True,
):
    """
    Obtain a bias field estimate in calibration space using
    the HCP's sebased approach.

    Parameters
    ----------
    calib_name: pathlib.Path
        Path to the calibration image.
    struct2asl: pathlib.Path
        Path to the .mat file which will be used to register
        the T1-estimated bias field to calibration image space.
    wmseg_name: pathlib.Path
        Path to a white matter segmentation for use in BBR
        registration.
    results_dir: pathlib.Path
        Directory in which to store results.
    fmapmag: pathlib.Path
        Path to the fieldmap magnitude image.
    fmapmagbrain: pathlib.Path
        Path to the brain-extracted fieldmap magnitude image.
    interpolation: int, default=3
        Order of interpolation to be used when registering
        the T1 image's bias field to the calibration image.
    force_refresh: bool, default=True
        Whether to recreate intermediate files if they already
        exist.

    Returns
    -------
    bias_field: Nifti1Image
    """
    # register fieldmap to t1 image
    fmap_reg_dir = results_dir / "fmap_registration"
    fmap_reg_dir.mkdir(exist_ok=True, parents=True)
    fmapmag_calib_name = fmap_reg_dir / "fmapmag_calibspc.nii.gz"
    struct2asl_reg = rt.Registration.from_flirt(struct2asl, src=t1_name, ref=calib_name)
    if not fmapmag_calib_name.exists() or force_refresh:
        register_fmap(
            fmapmag, fmapmagbrain, t1_name, t1_brain_name, fmap_reg_dir, wmseg_name
        )
        fmap_bbr = rt.Registration.from_flirt(
            fmap_reg_dir / "fmapmag2struct_bbr.mat",
            src=fmapmag,
            ref=t1_name,
        )
        fmap2calib = rt.chain(fmap_bbr, struct2asl_reg)
        fmap_calib = fmap2calib.apply_to_image(
            src=fmapmag, ref=calib_name, order=interpolation
        )
        nb.save(fmap_calib, fmapmag_calib_name)
    # get gray matter mask for sebased
    gm_seg_name = results_dir / "gm_mask.nii.gz"
    if not gm_seg_name.exists() or force_refresh:
        gm_seg = generate_tissue_mask_in_ref_space(
            aparc_aseg, calib_name, "gm", struct2asl, order=0
        )
        nb.save(gm_seg, gm_seg_name)
    # get brain mask
    brain_mask = results_dir / "brain_mask.nii.gz"
    if not brain_mask.exists() or force_refresh:
        t1_brain_img = nb.load(t1_brain_name)
        t1_mask = nb.nifti1.Nifti1Image(
            np.where(t1_brain_img.get_fdata() > 0, 1.0, 0.0), affine=t1_brain_img.affine
        )
        aslt1_mask = struct2asl_reg.apply_to_image(t1_mask, calib_name, order=0)
        aslt1_mask = nb.nifti1.Nifti1Image(
            np.where(aslt1_mask.get_fdata() > 0.5, 1.0, 0.0), affine=aslt1_mask.affine
        )
        nb.save(aslt1_mask, brain_mask)
    # get sebased bias
    bias_name = results_dir / "sebased_bias_dil.nii.gz"
    if not bias_name.exists() or force_refresh:
        sebased_cmd = [
            "get_sebased_bias_asl",
            "-i",
            calib_name,
            "-f",
            fmapmag_calib_name,
            "-m",
            brain_mask,
            "-o",
            results_dir,
            "--tissue_mask",
            gm_seg_name,
            "--debug",
        ]
        subprocess.run(sebased_cmd, check=True)
    dilall_name = results_dir / "sebased_bias_dilall.nii.gz"
    dilall_cmd = ["fslmaths", bias_name, "-dilall", dilall_name]
    subprocess.run(dilall_cmd, check=True)
    bias_field = nb.load(dilall_name)
    return bias_field


METHODS = ("calib", "t1", "sebased")
BIAS_ESTIMATION = {
    "calib": bias_estimation_calib,
    "t1": bias_estimation_t1,
    "sebased": bias_estimation_sebased,
}


def bias_estimation(calib_name, method, **kwargs):
    """
    Wrapper for the different bias field estimation
    functions.

    Parameters
    ----------
    calib_name: pathlib.Path
    method: str

    Returns
    -------
    bias_field: Nifti1Image

    There are also several kwargs which will be passed on
    to the wrapped functions depending on the chosen bias
    estimation method. These are documented briefly below.
    kwargs
    ------
    method == "calib"
        n/a
    method == "t1"
        fslanatdir: pathlib.Path
        struct2asl: pathlib.Path
        interpolation: int
    method == "sebased"
        fslanatdir: pathlib.Path
        struct2asl: pathlib.Path
        interpolation: int
        wmseg_name: pathlib.Path
        results_dir: pathlib.Path
        fmapmag: pathlib.Path
        fmapmagbrain: pathlib.Path
    """
    # check the calibration image exists
    calib_name = Path(calib_name).resolve(strict=True)
    # assert the bias estimation method is one of those supported
    assert method in METHODS, f"{method} is not a supported method."
    # call appropriate bias estimation function
    try:
        return BIAS_ESTIMATION[method](calib_name, **kwargs)
    except TypeError as e:
        raise TypeError(
            f"Incorrect arguments for {method} " f"based bias estimation.\n{e}"
        )
