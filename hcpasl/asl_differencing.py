"""
Functions to perform the tag-control differencing of multi-band 
ASL data which has been rescaled and motion corrected.

The rescaling is accounted for by using a GLM to perform the 
tag-control subtraction, similar to the method described in 
'A framework for motion correction of background suppressed
arterial spin labeling perfusion images acquired with
simultaneous multi‐slice EPI', Y. Suzuki, T.W. Okell, M.A. 
Chappell, M.J.P. van Osch
"""

import nibabel as nb
import numpy as np


def tag_control_differencing(series, scaling_factors, betas_dir, mask=None):
    """
    Perform tag-control differencing of a scaled ASL sequence.

    The method follows that in [1]_.

    Parameters
    ----------
    series : pathlib.Path
        Path to the ASL series we wish to difference.
    scaling_factors : pathlib.Path
        Path to the scaling factors used to de-band the ASL series.
    betas_dir : pathlib.Path
        Path to the directory in which to save the results.
    mask : pathlib.Path:
        operate only within mask

    .. [1] Suzuki, Yuriko, et al. "A framework for motion
       correction of background suppressed arterial spin labeling
       perfusion images acquired with simultaneous multi‐slice
       EPI." Magnetic resonance in medicine 81.3 (2019):
       1553-1565.
    """
    # create output directory
    betas_dir.mkdir(exist_ok=True, parents=True)

    # load motion- and distortion- corrected data, Y_moco
    Y_moco = nb.load(series)

    # load registered scaling factors, S_st
    S_st = nb.load(scaling_factors)

    # calculate X_perf = X_tc * S_st
    X_tc = np.ones((1, 1, 1, 86)) * 0.5
    X_tc[0, 0, 0, 0::2] = -0.5
    X_perf = X_tc * S_st.get_fdata()

    # split X_perf and Y_moco into even and odd indices
    X_odd = X_perf[:, :, :, 1::2]
    X_even = X_perf[:, :, :, 0::2]
    Y_odd = Y_moco.get_fdata()[:, :, :, 1::2]
    Y_even = Y_moco.get_fdata()[:, :, :, 0::2]

    # ignore voxels where below would lead to dividing by zero
    nonzero_mask = np.abs(X_odd - X_even) > 1e-6
    mask_name = betas_dir / "difference_mask.nii.gz"
    nb.save(nb.Nifti1Image(nonzero_mask.astype(int), affine=Y_moco.affine), mask_name)

    # only perform calculation within the provided mask
    if mask is not None:
        mask = nb.load(mask).get_fdata()
        if mask.ndim == 3:
            mask = mask[..., np.newaxis]
        nonzero_mask = np.logical_and(nonzero_mask, mask)
        mask_name = betas_dir / "combined_mask.nii.gz"
        nb.save(
            nb.Nifti1Image(nonzero_mask.astype(int), affine=Y_moco.affine), mask_name
        )

    # calculate B_perf and B_baseline at voxels within mask
    B_perf, B_baseline = np.zeros_like(X_odd), np.zeros_like(X_even)
    B_perf[nonzero_mask] = (Y_odd[nonzero_mask] - Y_even[nonzero_mask]) / (
        X_odd[nonzero_mask] - X_even[nonzero_mask]
    )
    B_baseline[nonzero_mask] = (
        X_odd[nonzero_mask] * Y_even[nonzero_mask]
        - X_even[nonzero_mask] * Y_odd[nonzero_mask]
    ) / (X_odd[nonzero_mask] - X_even[nonzero_mask])

    # save both images
    B_perf_name = betas_dir / "beta_perf.nii.gz"
    B_perf_img = nb.nifti1.Nifti1Image(B_perf, affine=Y_moco.affine)
    nb.save(B_perf_img, B_perf_name)
    B_baseline_name = betas_dir / "beta_baseline.nii.gz"
    B_baseline_img = nb.nifti1.Nifti1Image(B_baseline, affine=Y_moco.affine)
    nb.save(B_baseline_img, B_baseline_name)
