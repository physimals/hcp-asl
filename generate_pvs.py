import os.path as op 
import os

import numpy as np 
import nibabel as nib

import toblerone as tob
from toblerone.pvestimation import estimators
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

def extract_fs_pvs(aparcseg, ref): 

    superfactor = 2
    ref_spc = rt.ImageSpace(ref)
    high_spc = ref_spc.resize_voxels(1/superfactor, 'ceil')

    src_spc = nib.load(aparcseg)
    aparcseg = src_spc.get_fdata().astype(np.int32)
    src_spc = rt.ImageSpace(src_spc)

    # Need cortex_*, FAST_*, and individual structures 
    pv_dict = {}
    cortex = np.zeros((np.prod(aparcseg.shape), 3), dtype=np.int8)

    for label in np.unique(aparcseg):
        tissue = FS_LUT.get(label)
        if tissue: 
            mask = (aparcseg == label) 
            if tissue == "WM":
                cortex[mask.flatten(),1] = 1
            elif tissue == "GM":
                cortex[mask.flatten(),0] = 1 
            elif tissue == "CSF":
                pass 
            else: 
                pv_dict[tissue] = mask.astype(np.int8)
        else: 
            print("Did not assign", label)

    reg = rt.Registration.identity()
    cortex = src_spc.make_nifti(cortex.reshape(*aparcseg.shape, 3))
    cortex = reg.apply_to_image(cortex, high_spc, order=1)
    cortex = _sumArrayBlocks(cortex, 3*[superfactor] + [1]) / (superfactor**3)
    cortex = cortex.reshape(-1,3)
    cortex[:,2] = np.maximum(0, 1 - cortex[:,:2].sum(1))
    cortex[cortex[:,2] < 1e-2, 2] = 0 
    cortex /= cortex.sum(1)[:,None]
    cortex = cortex.reshape(*ref_spc.size, 3)

    keys, values = list(pv_dict.keys()), list(pv_dict.values())
    subcorts = src_spc.make_nifti(np.stack(values, axis=-1))
    subcorts = reg.apply_to_image(subcorts, high_spc, order=1)
    subcorts = _sumArrayBlocks(subcorts, 3*[superfactor] + [1]) / (superfactor**3)
    subcorts = np.moveaxis(subcorts, 3, 0)
    pv_dict = dict(zip(keys, subcorts))

    pv_dict['cortex_GM'] = cortex[...,0]
    pv_dict['cortex_WM'] = cortex[...,1]
    pv_dict['cortex_nonbrain'] = cortex[...,2]
    pv_dict['FAST_GM'] = pv_dict['cortex_GM']
    pv_dict['FAST_WM'] = pv_dict['cortex_WM']
    pv_dict['FAST_CSF'] = pv_dict['cortex_nonbrain']

    result = estimators.stack_images(pv_dict)

    # Squash small values (looks like dodgy output otherwise, but
    # is actually just an artefact of resampling)
    pvmin = result.min()
    if pvmin > 1e-9:
        result[pvs < 1e-6] = 0 

    ref_spc.save_image(result, "testdata/stacked.nii.gz")


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

if __name__ == "__main__":


    fsdir = 'testdata/fs'
    aparcseg = op.join(fsdir, 'mri/aseg.mgz')
    struct2ref = np.eye(4)
    ref = op.join('testdata/T1low.nii.gz')
    extract_fs_pvs(aparcseg, ref)