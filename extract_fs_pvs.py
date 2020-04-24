#! /usr/bin/env python3
"""Script to extract PVs from FS binary volumetric segmentation"""

import pathlib
import sys 
import argparse
import multiprocessing as mp 

import numpy as np 
import nibabel as nib

import regtools as rt 

FS_LUT = {
    # Left hemisphere
    2 : "WM",
    3 : "GM",
    4 : "CSF",
    5 : "CSF",
    26 : "L_Accu",
    18 : "L_Amyg",
    11 : "L_Caud",
    17 : "L_Hipp",
    13 : "L_Pall",
    12 : "L_Puta",
    9 : "L_Thal",
    10 : "L_Thal",

    24 : "CSF",

    # Right hemisphere 
    41 : "WM",
    42 : "GM",
    43 : "CSF",
    44 : "CSF",
    58 : "R_Accu",
    54 : "R_Amyg",
    50 : "R_Caud",
    53 : "R_Hipp",
    52 : "R_Pall",
    51 : "R_Puta",
    48 : "R_Thal",
    49 : "R_Thal"
}

def extract_fs_pvs(aseg, t1, asl, superfactor, cores): 

    asl_spc = rt.ImageSpace(asl)
    t1_spc = rt.ImageSpace(t1)
    ref_spc = t1_spc.resize_voxels(asl_spc.vox_size / t1_spc.vox_size)
    high_spc = ref_spc.resize_voxels(1/superfactor, 'ceil')

    src_spc = nib.load(aseg)
    aseg = src_spc.get_fdata().astype(np.int32)
    src_spc = rt.ImageSpace(src_spc)

    # Need cortex_*, FAST_*, and individual structures 
    pv_dict = {}
    cortex = np.zeros((np.prod(aseg.shape), 3), dtype=np.float32)

    for label in np.unique(aseg):
        tissue = FS_LUT.get(label)
        if tissue: 
            mask = (aseg == label) 
            if tissue == "WM":
                cortex[mask.flatten(),1] = 1
            elif tissue == "GM":
                cortex[mask.flatten(),0] = 1 
            elif tissue == "CSF":
                pass 
            else: 
                pv_dict[tissue] = mask.astype(np.float32)
        else: 
            print("Did not assign", label)

    reg = rt.Registration.identity()
    cortex = src_spc.make_nifti(cortex.reshape(*aseg.shape, 3))
    cortex = reg.apply_to_image(cortex, high_spc, order=1, cores=cores)
    cortex = _sumArrayBlocks(cortex, 3*[superfactor] + [1]) / (superfactor**3)
    cortex = cortex.reshape(-1,3)
    cortex[:,2] = np.maximum(0, 1 - cortex[:,:2].sum(1))
    cortex[cortex[:,2] < 1e-2, 2] = 0 
    cortex /= cortex.sum(1)[:,None]
    cortex = cortex.reshape(*ref_spc.size, 3)

    keys, values = list(pv_dict.keys()), list(pv_dict.values())
    subcorts = src_spc.make_nifti(np.stack(values, axis=-1))
    del pv_dict
    subcorts = reg.apply_to_image(subcorts, high_spc, order=1, cores=cores)
    subcorts = _sumArrayBlocks(subcorts, 3*[superfactor] + [1]) / (superfactor**3)
    subcorts = np.moveaxis(subcorts, 3, 0)
    pv_dict = dict(zip(keys, subcorts))

    pv_dict['cortex_GM'] = cortex[...,0]
    pv_dict['cortex_WM'] = cortex[...,1]
    pv_dict['cortex_nonbrain'] = cortex[...,2]
    result = stack_images(pv_dict)

    # Squash small values (looks like dodgy output otherwise, but
    # is actually just an artefact of resampling)
    pvmin = result.min()
    if pvmin > 1e-9:
        result[result < 1e-6] = 0 

    return ref_spc.make_nifti(result.reshape((*ref_spc.size, 3))) 


def _sumArrayBlocks(array, factor):
    """Sum sub-arrays of a larger array, each of which is sized according to factor. 
    The array is split into smaller subarrays of size given by factor, each of which 
    is summed, and the results returned in a new array, shrunk accordingly. 

    Args:
        array: n-dimensional array of data to sum
        factor: n-length tuple, size of sub-arrays to sum over

    Returns:
        array of size array.shape/factor, each element containing the sum of the 
            corresponding subarray in the input
    """

    if len(factor) != len(array.shape):
        raise RuntimeError("factor must be of same length as number of dimensions")

    if np.any(np.mod(factor, np.ones_like(factor))):
        raise RuntimeError("factor must be of integer values only")

    factor = [ int(f) for f in factor ]

    outshape = [ int(s/f) for (s,f) in zip(array.shape, factor) ]
    out = np.copy(array)

    for dim in range(3):
        newshape = [0] * 4

        for d in range(3):
            if d < dim: 
                newshape[d] = outshape[d]
            elif d == dim: 
                newshape[d+1] = factor[d]
                newshape[d] = outshape[d]
            else: 
                newshape[d+1] = array.shape[d]

        newshape = newshape + list(array.shape[3:])
        out = np.sum(out.reshape(newshape), axis=dim+1)

    return out 


def stack_images(images):
    """
    Combine the results of estimate_all() into overall PV maps
    for each tissue. Note that the below logic is entirely specific 
    to the surfaces produced by FreeSurfer, FIRST and how they may be
    combined with FAST estimates. If you're going off-piste anywhere else
    then you probably DON'T want to re-use this logic. 

    Args: 
        dictionary of PV maps, keyed as follows: all FIRST subcortical
        structures named by their FIST convention (eg L_Caud); FAST's estimates
        named as FAST_CSF/WM/GM; cortex estimates as cortex_GM/WM/non_brain

    Returns: 
        single 4D array of PVs, arranged GM/WM/non-brain in the 4th dim
    """

    # The logic is as follows: 
    # Initialise everything as non-brain
    # Write in PV estimates from the cortical surfaces. This sets the cortical GM
    # and also all subcortical tissue as WM 
    # Layer in subcortical GM from each individual FIRST structure (except brain stem). 
    # After adding in each structure's GM, recalculate CSF as either the existing amount, 
    # or reduce the CSF estimate if the total (GM + CSF) > 1
    # If there is any remainder unassigned in the voxel, set that as WM

    # To summarise, the tissues are stacked up as: 
    # All CSF 
    # Cortical GM
    # Then within the subcortex only: 
        # All subcortical volume set as WM 
        # Subcortical CSF fixed using FAST 
        # Add in subcortical GM for each structure, reducing CSF if required 
        # Set the remainder 1 - (GM+CSF) as WM in voxels that were updated 

    # Pop the cortex estimates and initialise output as all CSF
    ctxgm = images.pop('cortex_GM').flatten()
    ctxwm = images.pop('cortex_WM').flatten()
    ctxnon = images.pop('cortex_nonbrain').flatten()
    shape = (*ctxgm.shape, 3)
    ctx = np.vstack((ctxgm, ctxwm, ctxnon)).T
    out = np.zeros_like(ctx)
    out[:,2] = 1

    # Then write in cortex estimates from all voxels
    # that contain either WM or GM (on the ctx image)
    mask = np.logical_or(ctx[:,0], ctx[:,1])
    out[mask,:] = ctx[mask,:]

    # Sanity checks: total tissue PV in each vox should sum to 1
    # assert np.all(out[to_update,0] <= GM_threshold), 'Some update voxels have GM'
    assert np.all(np.abs(out.sum(1) - 1) < 1e-6), 'Voxel PVs do not sum to 1'

    # For each subcortical structure, create a mask of the voxels which it 
    # relates to. The following operations then apply only to those voxels 
    # All subcortical structures interpreted as pure GM 
    # Update CSF to ensure that GM + CSF in those voxels < 1 
    # Finally, set WM as the remainder in those voxels.
    for s in images.values():
        smask = (s.flatten() > 0)
        out[smask,0] = np.minimum(1, out[smask,0] + s.flatten()[smask])
        out[smask,2] = np.minimum(out[smask,2], 1 - out[smask,0])
        out[smask,1] = np.maximum(1 - (out[smask,0] + out[smask,2]), 0)

    # Final sanity check, then rescaling so all voxels sum to unity. 
    assert (out > -1e-6).all()
    out[out < 0] = 0 
    sums = out.sum(1)
    assert (np.abs(sums - 1) < 1e-6).all()
    out = out / sums[:,None]

    return out.reshape(shape)


if __name__ == "__main__":

    usage = """
    usage: generate_pvs --fsdir <dir> --t1 <path> --asl <path> --out <path> [--stack]

    Extracts PV estimates from FS' binary volumetric segmentation. Saves output 
    in T1 space, with voxels resized to match ASL resolution. 
    
    --aseg: path to aseg.mgz in FS output directory
    --t1: path to T1 image for alignment
    --asl: path to ASL image
    --out: path to save output at (with suffix _GM, _WM, _CSF, unless --stack)
    --stack: return 4D volume, stacked GM,WM,CSF in last dimension (default false)
    --super: super-sampling level for intermediate resampling (default 4)
    --cores: CPU cores to use (default max)
    """

    # cmd = "--aseg testdata/fs/mri/aseg.mgz --t1 testdata/T1.nii.gz --asl testdata/T1low.nii.gz --out testdata/stacked.nii.gz"
    # sys.argv = cmd.split()

    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument("--aseg", required=True)
    parser.add_argument("--t1", required=True)
    parser.add_argument("--asl", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--stack", action="store_true")
    parser.add_argument("--super", default=4, type=int)
    parser.add_argument("--cores", default=mp.cpu_count(), type=int)

    args = parser.parse_args(sys.argv)
    pvs = extract_fs_pvs(args.aseg, args.t1, args.asl, args.super, args.cores)

    if args.stack: 
        nib.save(pvs, args.out)

    else: 
        opath = pathlib.Path(args.out)
        spc = rt.ImageSpace(pvs)
        pvs = pvs.dataobj
        for idx,tiss in enumerate(["GM", "WM", "CSF"]):
            n = opath.parent.joinpath(opath.stem.rsplit('.', 1)[0] + f"_{tiss}" + "".join(opath.suffixes))
            spc.save_image(pvs[...,idx], n.as_posix())