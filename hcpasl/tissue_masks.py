import regtricks as rt
import nibabel as nb

import numpy as np
import scipy

TISSUE_LABELS = {
    "wm": (2, 41),
    "allwm":(2, 41, 77, *list(range(251, 256)), 7, 46),
    "csf": (4, 43),
    "allvent": (4, 43, 5, 14, 44, 15, 72, 31, 63, 0, 24)
}

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
    with open(LUT_name) as f:
        labels = [int(l.split(" ")[0]) for l in f.readlines()[1::2]]
    return labels

def generate_tissue_mask(aparc_aseg, tissue, erode=False):
    """
    Generate a tissue mask from FreeSurfer's aparc+aseg.

    Parameters
    ----------
    aparc_aseg: Path to FS's T1-space aparac_aseg output.
    tissue: str, tissue of interest.
    erode: bool, whether to erode the initial mask or not, 
            default is False.

    Returns
    -------
    mask: np.array logical mask of tissue of interest
    """
    # load aparc_aseg
    aseg = nb.load(aparc_aseg).get_fdata()

    # create mask
    mask = np.zeros_like(aseg)

    # get appropriate labels
    if tissue == "gm":
        labels = (*TISSUE_LABELS["allwm"], *TISSUE_LABELS["allvent"])
    else:
        labels = TISSUE_LABELS[tissue]

    # iterate through labels
    for label in labels:
        mask = np.where(aseg==label, 1., mask)
    
    # invert mask if tissue==gm
    if tissue == "gm":
        mask = np.where(mask==0., 1., 0.)

    # potential round of eroding
    if erode:
        mask = scipy.ndimage.morphology.binary_erosion(mask)
    
    return mask

def generate_tissue_mask_in_ref_space(aparc_aseg, 
                                      ref_img,
                                      tissue, 
                                      struct2ref=None, 
                                      superfactor=True, 
                                      order=3, 
                                      threshold=0.8,
                                      erode=False):
    """
    Generate a tissue mask from FreeSurfer's aparc+aseg in the 
    space of the given reference image.

    Parameters
    ----------
    tissue: Which tissue mask we would like, (wm, allwm, csf, 
            allvent, gm).
    aparc_aseg: Path to FS's T1-space aparc_aseg output.
    ref_img: Path to the reference space image.
    struct2ref: Path to the registration from structural to 
                the reference space, default is None (assumes 
                the reference image is already in the desired 
                space).
    order: Order of interpolation to be used when applying 
           registration, default is 3.
    threshold: Threshold to use to re-binarise the tissue 
               segmentation after registration, default is 0.8.
    erode: bool, whether to erode the initial mask (in original 
            structural space) or not, default is False.
    
    Returns
    -------
    mask: Nifti1Image
    """

    # load reference image
    ref_img = nb.load(ref_img)

    # get T1 space tissue mask
    mask = generate_tissue_mask(aparc_aseg, tissue, erode=erode)

    # resample to reference image
    if struct2ref:
        reg = rt.Registration.from_flirt(struct2ref,
                                         src=aparc_aseg,
                                         ref=ref_img)
    else:
        reg = rt.Registration.identity()
    mask = reg.apply_to_array(
        data=mask, src=aparc_aseg, ref=ref_img, 
        superfactor=superfactor, order=order
    )

    # re-binarise
    mask = np.where(mask>=threshold, 1., 0.)

    # create nifti
    mask = nb.nifti1.Nifti1Image(mask, affine=ref_img.affine)
    return mask
