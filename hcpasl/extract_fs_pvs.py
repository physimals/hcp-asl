#! /usr/bin/env python3
"""Script to extract PVs from FS binary volumetric segmentation"""

import argparse
import pathlib

import nibabel as nib
import numpy as np
import regtricks as rt
from toblerone.pvestimation import cortex as estimate_cortex

# Map from FS aseg labels to tissue types
# NB cortex is ignored
FS_LUT = {
    # Left hemisphere
    2: "WM",
    3: "GM",
    28: "WM",  # Left ventral DC
    30: "WM",  # Left vessel
    31: "GM",  # Left choroid plexus
    26: "GM",  # L_Accu
    18: "GM",  # L_Amyg
    11: "GM",  # L_Caud
    17: "GM",  # L_Hipp
    13: "GM",  # L_Pall
    12: "GM",  # L_Puta
    9: "GM",  # L_Thal
    10: "GM",  # L_Thal
    77: "WM",
    78: "WM",
    79: "WM",
    # Right hemisphere
    41: "WM",
    42: "GM",
    60: "WM",  # Right ventral DC
    62: "WM",  # Right vessel
    63: "GM",  # Right choroid plexus
    58: "GM",  # R_Accu
    54: "GM",  # R_Amyg
    50: "GM",  # R_Caud
    53: "GM",  # R_Hipp
    52: "GM",  # R_Pall
    51: "GM",  # R_Puta
    48: "GM",  # R_Thal
    49: "GM",  # R_Thal
    # Left cerebellum
    7: "WM",
    8: "GM",
    # Right cerebellum
    46: "WM",
    47: "GM",
    16: "WM",  # Brainstem
    85: "WM",  # Optic chiasm
}


def extract_fs_pvs(aparcseg, surf_dict, ref_spc, ref2struct=None, cores=1):
    """
    Extract and layer PVs according to tissue type, taken from a FS aparc+aseg.
    Results are stored in ASL-gridded T1 space.
    Args:
        aparcseg: path to aparc+aseg file
        surf_dict: dict with LWS/LPS/RWS/RPS keys, paths to those surfaces
        ref_spc: space in which to estimate (ie, ASL-gridded T1)
        ref2struct: FLIRT registration between reference and T1w image 
        cores: number CPU cores to use, default is 1
    Returns:
        nibabel Nifti object
    """

    ref_spc = rt.ImageSpace(ref_spc)
    aseg_spc = nib.load(aparcseg)
    aseg = aseg_spc.get_fdata().astype(int)
    aseg_spc = rt.ImageSpace(aseg_spc)

    # If not provided, an identity transform is used
    if ref2struct:
        struct2ref_reg = rt.Registration.from_flirt(
            ref2struct, ref_spc, aseg_spc
        ).inverse()
    else:
        struct2ref_reg = rt.Registration.identity()

    # Extract PVs from aparcseg segmentation, ignoring cortex.
    non_cortex_pvs = np.zeros((*aseg_spc.size, 2), dtype=np.float32)
    for k, t in FS_LUT.items():
        m = aseg == k
        idx = ["GM", "WM"].index(t)
        non_cortex_pvs[m, idx] = 1

    # Super-resolution resampling for the vol_pvs, a la applywarp.
    # 0: GM, 1: WM, always in the LAST dimension of an array
    non_cortex_pvs = struct2ref_reg.apply_to_array(
        non_cortex_pvs, aseg_spc, ref_spc, order=1, superfactor=5, cores=cores
    )

    ref_spc.save_image(non_cortex_pvs, "non_cortex.nii.gz")

    # Estimate cortical PVs
    cortex = estimate_cortex(
        ref=ref_spc,
        struct2ref=struct2ref_reg.src2ref,
        cores=cores,
        **surf_dict,
    )

    out = non_cortex_pvs[..., :2].copy()

    # Voxels where toblerone has identified the cortex
    ctx = cortex[..., 0] > 0.01

    # total brain PV in these voxels (GM or WM)
    ctx_brain_pv = np.maximum(cortex[ctx, :2].sum(-1), non_cortex_pvs[ctx, :2].sum(-1))

    # Layer in cortical GM (sum on top of existing GM)
    out[ctx, 0] = np.minimum(cortex[ctx, 0] + out[ctx, 0], 1)

    # In those voxels, the total brain PV be as predicted by Toblerone,
    # so update the WM value based on this
    out[ctx, 1] = np.maximum(ctx_brain_pv - out[ctx, 0], 0)

    # Sanity checks
    assert (out >= 0).all(), "Negative PV"
    assert (out <= 1).all(), "PV > 1"
    assert out.sum(-1).max() <= 1.001, "PV sum > 1"

    ref_spc.save_image(out, "asl_pvs.nii.gz")
    ref_spc.save_image(out.sum(-1), "gmwm.nii.gz")

    return ref_spc.make_nifti(out)


def main():
    desc = """
    Generate PV estimates for a reference voxel grid. FS' volumetric 
    segmentation is used for the subcortex and a surface-based method
    is used for the cortex (toblerone). 
    """

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument(
        "--aparcseg",
        required=True,
        help="path to volumetric FS segmentation (aparc+aseg.mgz)",
    )
    parser.add_argument("--LWS", required=True, help="path to left white surface")
    parser.add_argument("--LPS", required=True, help="path to left pial surface")
    parser.add_argument("--RPS", required=True, help="path to right pial surface")
    parser.add_argument("--RWS", required=True, help="path to right white surface")
    parser.add_argument(
        "--ref", required=True, help="path to image defining reference grid for PVs"
    )
    parser.add_argument("--out", required=True, help="path to save output")
    parser.add_argument(
        "--stack", action="store_true", help="stack output into single 4D volume"
    )
    parser.add_argument("--cores", default=1, type=int, help="CPU cores to use")

    args = parser.parse_args()
    surf_dict = dict([(k, getattr(args, k)) for k in ["LWS", "LPS", "RPS", "RWS"]])

    pvs = extract_fs_pvs(
        aparcseg=args.aparcseg, surf_dict=surf_dict, ref_spc=args.ref, cores=args.cores
    )

    if args.stack:
        nib.save(pvs, args.out)

    else:
        opath = pathlib.Path(args.out)
        spc = rt.ImageSpace(pvs)
        pvs = pvs.dataobj
        for idx, tiss in enumerate(["GM", "WM", "CSF"]):
            n = opath.parent.joinpath(
                opath.stem.rsplit(".", 1)[0] + f"_{tiss}" + "".join(opath.suffixes)
            )
            spc.save_image(pvs[..., idx], n.as_posix())
