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

from .utils import load_json, update_json
import nibabel as nb
from pathlib import Path
import subprocess
import numpy as np

def tag_control_differencing(series, scaling_factors, betas_dir, subject_dir, outdir="hcp_asl"):
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
    subject_dir : pathlib.Path
        Path to the subject's base directory.
    outdir : str
        Name of the main results directory. Default is 'hcp_asl'.

    .. [1] Suzuki, Yuriko, et al. "A framework for motion 
       correction of background suppressed arterial spin labeling 
       perfusion images acquired with simultaneous multi‐slice 
       EPI." Magnetic resonance in medicine 81.3 (2019): 
       1553-1565.
    """
    # load subject's json
    json_dict = load_json(subject_dir/outdir)

    # load motion- and distortion- corrected data, Y_moco
    Y_moco = nb.load(series)

    # load registered scaling factors, S_st
    S_st = nb.load(scaling_factors)

    # calculate X_perf = X_tc * S_st
    X_tc = np.ones((1, 1, 1, 86)) * 0.5
    X_tc[0, 0, 0, 0::2] =  -0.5
    X_perf = X_tc * S_st.get_fdata()

    # split X_perf and Y_moco into even and odd indices
    X_odd = X_perf[:, :, :, 1::2]
    X_even = X_perf[:, :, :, 0::2]
    Y_odd = Y_moco.get_data()[:, :, :, 1::2]
    Y_even = Y_moco.get_data()[:, :, :, 0::2]

    # calculate B_perf and B_baseline
    B_perf = (Y_odd - Y_even) / (X_odd - X_even)
    B_baseline = (X_odd*Y_even - X_even*Y_odd) / (X_odd - X_even)

    # save both images
    betas_dir.mkdir(exist_ok=True)
    B_perf_name = betas_dir / 'beta_perf.nii.gz'
    B_perf_img = nb.nifti1.Nifti1Image(B_perf,
                                       affine=Y_moco.affine)
    nb.save(B_perf_img, B_perf_name)
    B_baseline_name = betas_dir / 'beta_baseline.nii.gz'
    B_baseline_img = nb.nifti1.Nifti1Image(B_baseline,
                                           affine=Y_moco.affine)
    nb.save(B_baseline_img, B_baseline_name)

    # add B_perf_name to the json as will be needed in oxford_asl
    important_names = {
        'beta_perf': str(B_perf_name)
    }
    update_json(important_names, json_dict)