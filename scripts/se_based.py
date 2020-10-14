import subprocess
import argparse
from pathlib import Path
import shutil
import regtricks as rt
import nibabel as nb
import numpy as np
import scipy
import pandas as pd
from fsl.data.image import Image

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
    lut = pd.read_csv(LUT_name, header=None)
    labels = [int(row[0].split(" ")[0]) for _, row in lut[1::2].iterrows()]
    return labels

def se_based_bias_estimation():
    """
    This script seeks to replicate the SE-based bias estimation 
    in the HCP's ComputeSpinEchoBiasField.sh.

    Some modifications need to be made for the ASL pipeline so, 
    for now, we shall use the script below. Necessary 
    modifications include the use of an already pre-computed 
    SpinEchoMean.nii.gz (from topup), the re-sampling of 
    ribbon.mgz and wmparc.mgz into our ASL-gridded T1 space, 
    and the use of our proton density weighted images rather 
    than the GRE.nii.gz in the original script.
    """
    # argument handling
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        help="Image from which we wish to estimate the bias field."
    )
    parser.add_argument(
        '--asl',
        help="ASL series to which we wish to apply the bias field."
    )
    parser.add_argument(
        "-f",
        "--fmapmag",
        help="Fieldmap magnitude image from topup."
    )
    parser.add_argument(
        "-m",
        "--mask",
        help="Brain mask."
    )
    parser.add_argument(
        '--wmparc',
        help="wmparc.mgz from FreeSurfer"
    )
    parser.add_argument(
        '--ribbon',
        help="ribbon.mgz from FreeSurfer"
    )
    parser.add_argument(
        "--cortical",
        help="Filename for FreeSurfer's Cortical Lable Table"
    )
    parser.add_argument(
        "--subcortical",
        help="Filename for FreeSurfer's Subcortical Lable Table"
    )
    parser.add_argument(
        '-o',
        '--outdir',
        help="Output directory for results."
    )
    parser.add_argument(
        '--debug',
        help="If this argument is specified, all intermediate files "
            +"will be saved for inspection.",
        action='store_true'
    )
    parser.add_argument(
        '--wm_mask',
        help="Filename for white matter mask to use this instead of "
            +"a gray matter mask derived from FreeSurfer outputs.",
        default=None
    )

    args = parser.parse_args()

    m0_name = args.input
    asl_name = args.asl
    sem_name = args.fmapmag
    mask_name = args.mask
    wmparc_name = args.wmparc
    ribbon_name = args.ribbon
    cLUT_name = args.cortical
    scLUT_name = args.subcortical
    outdir = Path(args.outdir)
    debug = args.debug
    wm_mask = args.wm_mask

    # create output directory
    outdir.mkdir(exist_ok=True)

    # load images
    m0_img, sem_img, mask_img = [
        Image(name) for name in (m0_name, sem_name, mask_name)
    ]

    # find ratio between SpinEchoMean and M0
    SEdivM0 = np.where(m0_img.data!=0, (sem_img.data/m0_img.data), 0)
    if debug:
        SEdivM0_name = str(outdir/'SEdivM0.nii.gz')
        SEdivM0_img = Image(SEdivM0, header=m0_img.header)
        SEdivM0_img.save(SEdivM0_name)
        
    # apply mask to ratio
    SEdivM0_brain = SEdivM0 * mask_img.data
    if debug:
        SEdivM0_brain_name = str(outdir/'SEdivM0_brain.nii.gz')
        SEdivM0_brain_img = Image(SEdivM0_brain, header=m0_img.header)
        SEdivM0_brain_img.save(SEdivM0_brain_name)

    # get summary stats for thresholding
    nanned_temp = np.where(mask_img.data==0, np.nan, SEdivM0_brain)
    median, std = [np.nanmedian(nanned_temp), np.nanstd(nanned_temp)]
    if debug:
        print(np.array(median), np.array(std))
        savenames = [str(outdir/f'ratio_{stat}.txt') for stat in ('median', 'std')]
        [np.savetxt(name, [val]) for name, val in zip(savenames, (median, std))]

    # apply thresholding
    lower, upper = [median - std, median + std]
    print(lower, upper)
    SEdivM0_brain_thr = np.where(
        np.logical_and(SEdivM0_brain >= lower, SEdivM0_brain <= upper),
        SEdivM0_brain,
        0
    )
    if debug:
        SEdivM0_brain_thr_name = str(outdir/'SEdivM0_brain_thr.nii.gz')
        SEdivM0_brain_thr_img = Image(SEdivM0_brain_thr, header=m0_img.header)
        SEdivM0_brain_thr_img.save(SEdivM0_brain_thr_name)

    ### HCP pipeline does median dilation here but isn't used - skip for now ###

    # set sigma for smoothing used in HCPPipeline
    fwhm = 5
    sigma = fwhm / np.sqrt(8*np.log(2))
    # binarise and smooth the thresholded image
    SEdivM0_brain_thr_roi = np.where(SEdivM0_brain_thr>0, 1, 0).astype(np.float)
    SEdivM0_brain_thr_s5 = scipy.ndimage.gaussian_filter(SEdivM0_brain_thr, sigma=sigma)
    SEdivM0_brain_thr_roi_s5 = scipy.ndimage.gaussian_filter(SEdivM0_brain_thr_roi, sigma=sigma)
    SEdivM0_brain_bias = np.where(
        np.logical_and(SEdivM0_brain_thr_roi_s5!=0, mask_img.data==1),
        (SEdivM0_brain_thr_s5 / SEdivM0_brain_thr_roi_s5),
        0
    )
    if debug:
        savenames = [
                        str(outdir/f'SEdivM0_brain_{name}.nii.gz') 
                        for name in ('thr_roi', 'thr_s5', 'thr_roi_s5', 'bias')
                    ]
        images = [
            Image(array, header=m0_img.header)
            for array in (SEdivM0_brain_thr_roi, SEdivM0_brain_thr_s5, SEdivM0_brain_thr_roi_s5, SEdivM0_brain_bias)
        ]
        [image.save(savename) for image, savename in zip(images, savenames)]

    ### HCP pipeline does median dilation here but isn't used - skip for now ###

    # correct the SEFM image
    SpinEchoMean_brain_BC = np.where(
        np.logical_and(mask_img.data==1, SEdivM0_brain_bias!=0), 
        sem_img.data/SEdivM0_brain_bias, 
        0
    )
    if debug:
        SpinEchoMean_brain_BC_name = str(outdir/'SpinEchoMean_brain_BC.nii.gz')
        SpinEchoMean_brain_BC_img = Image(SpinEchoMean_brain_BC, header=m0_img.header)
        SpinEchoMean_brain_BC_img.save(SpinEchoMean_brain_BC_name)

    # get ratio between bias-corrected FM and M0 image
    SEBCdivM0_brain = np.where(
        np.logical_and(mask_img.data==1, SpinEchoMean_brain_BC!=0), 
        m0_img.data/SpinEchoMean_brain_BC, 
        0
    )
    if debug:
        SEBCdivM0_brain_name = str(outdir/'SEBCdivM0_brain.nii.gz')
        SEBCdivM0_brain_img = Image(SEBCdivM0_brain, header=m0_img.header)
        SEBCdivM0_brain_img.save(SEBCdivM0_brain_name)

    # find dropouts
    Dropouts = np.where(
        np.logical_and(SEBCdivM0_brain>0, SEBCdivM0_brain<0.6), 1, 0)
    Dropouts_inv = np.where(Dropouts==1, 0, 1)
    if debug:
        savenames = [str(outdir/f'{name}.nii.gz') for name in ('Dropouts', 'Dropouts_inv')]
        images = [Image(array, header=m0_img.header) for array in (Dropouts, Dropouts_inv)]
        [image.save(savename) for image, savename in zip(images, savenames)]
    
    if wm_mask:
        tissue_mask = rt.Registration.identity().apply_to_image(wm_mask, m0_name, order=0).get_fdata()
        if debug:
            savename = str(outdir/'WhiteMatter.nii.gz')
            image = Image(tissue_mask, header=m0_img.header)
            image.save(savename)
    else:
        # downsample wmparc and ribbon to ASL-gridded T1 resolution
        wmparc_aslt1, ribbon_aslt1 = [
            rt.Registration.identity().apply_to_image(name, m0_name, order=0)
            for name in (wmparc_name, ribbon_name)
        ]
        # parse LUTs
        c_labels, sc_labels = [parse_LUT(lut) for lut in (cLUT_name, scLUT_name)]
        cgm, scgm = [np.zeros(ribbon_aslt1.shape), np.zeros(wmparc_aslt1.shape)]
        cgm, scgm = [np.zeros(ribbon_aslt1.shape), np.zeros(wmparc_aslt1.shape)]
        for label in c_labels:
            cgm = np.where(ribbon_aslt1.get_fdata()==label, 1, cgm)
        for label in sc_labels:
            scgm = np.where(wmparc_aslt1.get_fdata()==label, 1, scgm)
        if debug:
            savenames = [
                str(outdir/f'{pre}GreyMatter.nii.gz') for pre in ('Cortical', 'Subcortical')
            ]
            images = [Image(array, header=m0_img.header) for array in (cgm, scgm)]
            [image.save(savename) for image, savename in zip(images, savenames)]
        
        # combine masks
        tissue_mask = np.where(np.logical_or(cgm==1, scgm==1), 1, 0)
        if debug:
            savename = str(outdir/'AllGreyMatter.nii.gz')
            image = Image(tissue_mask, header=m0_img.header)
            image.save(savename)
    
    # mask M0 image with both the tissue mask and Dropouts_inv mask
    M0_grey = np.where(np.logical_and(tissue_mask==1, Dropouts_inv==1), m0_img.data, 0).astype(np.float)
    M0_greyroi = np.where(M0_grey!=0, 1, 0).astype(np.float)
    M0_grey_s5, M0_greyroi_s5 = [
        scipy.ndimage.gaussian_filter(arr, sigma=sigma) for arr in (M0_grey, M0_greyroi)
    ]
    if debug:
        savenames = [str(outdir/f'M0_grey{part}.nii.gz') for part in ('', 'roi', '_s5', 'roi_s5')]
        images = [Image(array, header=m0_img.header) 
                  for array in (M0_grey, M0_greyroi, M0_grey_s5, M0_greyroi_s5)]
        [image.save(savename) for image, savename in zip(images, savenames)]

    # M0_bias_raw needs to undergo fslmaths' -dilall
    M0_bias_raw = np.where(
        np.logical_and(tissue_mask!=0, M0_greyroi_s5!=0), 
        M0_grey_s5/M0_greyroi_s5, 
        0
    )
    M0_bias_raw_name = str(outdir/'M0_bias_raw.nii.gz')
    M0_bias_raw_img = Image(M0_bias_raw, header=m0_img.header)
    M0_bias_raw_img.save(M0_bias_raw_name)
    dilall_cmd = ['fslmaths', M0_bias_raw_name, '-dilall', '-mas', mask_name, M0_bias_raw_name]
    subprocess.run(dilall_cmd, check=True)

    # reload dilalled-and-masked M0_bias_raw
    M0_bias_raw_img = Image(M0_bias_raw_name)
    M0_bias_raw = M0_bias_raw_img.data

    # refine bias field
    M0_bias_roi = np.where(M0_bias_raw>0, 1, 0).astype(np.float)
    M0_bias_raw_s5, M0_bias_roi_s5 = [
        scipy.ndimage.gaussian_filter(array, sigma=sigma) for array in (M0_bias_raw, M0_bias_roi)
    ]
    M0_bias = np.where(
        np.logical_and(mask_img.data!=0, M0_bias_roi_s5!=0), 
        M0_bias_raw_s5/M0_bias_roi_s5, 
        0
    )
    if debug:
        savenames = [
            str(outdir/f'M0_bias{part}.nii.gz')for part in ('roi', 'raw_s5', 'roi_s5', '')
        ]
        images = [
            Image(array, header=m0_img.header) 
            for array in (M0_bias_roi, M0_bias_raw_s5, M0_bias_roi_s5, M0_bias)
        ]
        [image.save(savename) for image, savename in zip(images, savenames)]
    
    # get summary stats
    nanned_temp = np.where(mask_img.data==0, np.nan, M0_bias)
    mean = np.nanmean(nanned_temp)

    # get sebased bias - should get ref also but leaving for now
    sebased_bias = M0_bias / mean
    sebased_bias_name = str(outdir/'sebased_bias.nii.gz')
    sebased_bias_img = Image(sebased_bias, header=m0_img.header)
    sebased_bias_img.save(sebased_bias_name)

    # apply 2 rounds of dilation to sebased_bias
    sebased_bias_dil_name = str(outdir/'sebased_bias_dil.nii.gz')
    bias_dil_cmd = ['fslmaths', sebased_bias_name, '-dilM', '-dilM', sebased_bias_dil_name]
    subprocess.run(bias_dil_cmd, check=True)

    # apply bias field to calibration and ASL images
    sebased_bias = Image(sebased_bias_dil_name)
    asl_img = Image(asl_name)
    calib_bc = np.where(sebased_bias.data!=0, m0_img.data/sebased_bias.data, m0_img.data)
    asl_bc = np.where(
        sebased_bias.data[..., np.newaxis]!=0, 
        asl_img.data/sebased_bias.data[..., np.newaxis], 
        asl_img.data
    )
    [
        Image(bc_arr, header=asl_img.header).save(str(outdir/f'{name}_secorr.nii.gz'))
        for name, bc_arr in zip(('calib0', 'tis'), (calib_bc, asl_bc))
    ]