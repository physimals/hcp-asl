import os
import subprocess as sp
import sys
import os.path as op 
import glob 
import tempfile 
from pathlib import Path
import multiprocessing as mp
import argparse

import regtricks as rt
import nibabel as nb
import scipy.ndimage
import numpy as np
from fsl.wrappers import fslmaths, bet
from fsl.data.image import Image
from scipy.ndimage import binary_fill_holes

from hcpasl.distortion_correction import (
    generate_gdc_warp, generate_topup_params, generate_fmaps, 
    generate_epidc_warp, register_fmap
)

def generate_asl2struct_initial(asl_vol0, struct, fsdir, reg_dir):
    """
    Generate the initial linear transformation between ASL-space and T1w-space
    using FS bbregister. This is further refined later on using asl_reg. Note
    that struct is required only for saving the output in the right convention, 
    it is not actually used by bbregister. 
    
    Args:
        asl_vol0: path to first volume of ASL 
        struct: path to T1w image (eg T1w_acdc_restore.nii.gz)
        fsdir: path to subject's FreeSurfer output directory 
        reg_dir: path to registration directory, for output 

    Returns: 
        n/a, file 'asl2struct_initial_bbr_fsl.mat' will be saved in reg_dir
    """

    # We need to do some hacky stuff to get bbregister to work...
    # Split the path to the FS directory into a fake $SUBJECTS_DIR
    # and subject_id. We temporarily set the environment variable
    # before making the call, and then revert back afterwards  
    new_sd, sid = op.split(fsdir)
    old_sd = os.environ.get('SUBJECTS_DIR')
    pwd = os.getcwd()
    orig_mgz = op.join(fsdir, 'mri', 'orig.mgz')

    # Run inside the regdir. Save the output in fsl format, by default 
    # this targets the orig.mgz, NOT THE T1 IMAGE ITSELF! 
    os.chdir(reg_dir)
    omat_path = op.join(reg_dir, "asl2orig_mgz_initial_bbr_fsl.mat")
    cmd = os.environ['SUBJECTS_DIR'] = new_sd
    cmd = f"$FREESURFER_HOME/bin/bbregister --s {sid} --mov {asl_vol0} --t1 "
    cmd += f"--reg asl2orig_mgz_initial_bbr.dat --fslmat {omat_path}"
    sp.run(cmd, shell=True)

    try:
        asl2orig_fsl = rt.Registration.from_flirt(omat_path, asl_vol0, orig_mgz)
    except RuntimeError as e:
        # final row != [0 0 0 1], round to 5 d.p. and try again
        print(e)
        print("Rounding to 5 d.p.")
        arr = np.loadtxt(omat_path)
        np.savetxt(omat_path, arr, fmt='%.5f')
        asl2orig_fsl = rt.Registration.from_flirt(omat_path, asl_vol0, orig_mgz)

    # Return to original working directory, and flip the FSL matrix to target
    # asl -> T1, not orig.mgz. Save output. 
    os.chdir(pwd)
    if old_sd:
        os.environ['SUBJECTS_DIR'] = old_sd
    asl2struct_fsl = asl2orig_fsl.to_flirt(asl_vol0, struct)
    np.savetxt(op.join(reg_dir, 'asl2struct_initial_bbr_fsl.mat'), asl2struct_fsl)

def generate_wmmask(aparc_aseg):
    """ 
    Generate binary WM mask in space of T1 image using FS aparc+aseg

    Args: 
        aparc_aseg: path to aparc_aseg in T1 space (not FS 256 1mm space!)

    Returns: 
        np.array logical WM mask in space of T1 image 
    """

    aseg_array = nb.load(aparc_aseg).get_data()
    wm = np.logical_or(aseg_array == 41, aseg_array == 2)
    return wm 

def generate_asl_mask(struct_brain, asl, asl2struct):
    """
    Generate brain mask in ASL space 

    Args: 
        struct_brain: path to T1 brain-extracted, ac_dc_restore_brain
        asl: path to ASL image 
        asl2struct: regtricks.Registration for asl to structural 

    Returns: 
        np.array, logical mask. 
    """

    brain_mask = (nb.load(struct_brain).get_data() > 0).astype(np.float32)
    asl_mask = asl2struct.inverse().apply_to_array(brain_mask, struct_brain, asl)
    asl_mask = binary_fill_holes(asl_mask > 0.25)
    return asl_mask

def binarise_image(image, threshold=0):
    """
    Binarise image above a threshold if given.

    Args:
        image: path to the image to be binarised
        threshold: voxels with a value below this will be zero and above will be one
    
    Returns:
        np.array, logical mask
    """
    image = Image(image)
    mask = (image.data>threshold).astype(np.float32)
    return mask

def create_ti_image(asl, tis, sliceband, slicedt, outname):
    """
    Create a 4D series of actual TIs at each voxel.

    Args:
        asl: path to image in the space we wish to create the TI series
        tis: list of TIs in the acquisition
        sliceband: number of slices per band in the acquisition
        slicedt: time taken to acquire each slice
        outname: path to which the ti image is saved
    
    Returns:
        n/a, file outname is created in output directory
    """

    asl_spc = rt.ImageSpace(asl)
    n_slice = asl_spc.size[2]
    slice_in_band = np.tile(np.arange(0, sliceband), 
                            n_slice//sliceband).reshape(1, 1, n_slice, 1)
    ti_array = np.array([np.tile(x, asl_spc.size) for x in tis]).transpose(1, 2, 3, 0)
    ti_array = ti_array + (slice_in_band * slicedt)
    rt.ImageSpace.save_like(asl, ti_array, outname)

def register_fmap(fmapmag, fmapmagbrain, s, sbet, out_dir, wm_tissseg):
    # create output directory
    out_dir = Path(out_dir)
    out_dir.mkdir(exist_ok=True)

    # get schedule
    fsldir = os.environ.get('FSLDIR')
    schedule = Path(fsldir)/'etc/flirtsch/bbr.sch'

    # set up commands
    init_xform = out_dir/'fmapmag2struct_init.mat'
    sec_xform = out_dir/'fmapmag2struct_sec.mat'
    bbr_xform = out_dir/'fmapmag2struct_bbr.mat'
    init_cmd = [
        'flirt',
        '-in', fmapmagbrain,
        '-ref', sbet,
        '-dof', '6',
        '-omat', init_xform
    ]
    sec_cmd = [
        'flirt',
        '-in', fmapmag,
        '-ref', s,
        '-dof', '6',
        '-init', init_xform,
        '-omat', sec_xform,
        '-nosearch'
    ]
    bbr_cmd = [
        'flirt',
        '-ref', s,
        '-in', fmapmag,
        '-dof', '6',
        '-cost', 'bbr',
        '-wmseg', wm_tissseg,
        '-init', sec_xform,
        '-omat', bbr_xform,
        '-schedule', schedule
    ]
    for cmd in (init_cmd, sec_cmd, bbr_cmd):
        sp.run(cmd, check=True)
    return str(bbr_xform)

def main():

    # argument handling
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--study_dir",
        help="Path of the base study directory.",
        required=True
    )
    parser.add_argument(
        "--sub_id",
        help="Subject number.",
        required=True
    )
    parser.add_argument(
        "-g",
        "--grads",
        help="Filename of the gradient coefficients for gradient"
            + "distortion correction.",
        required=True
    )
    parser.add_argument(
        "-t",
        "--target",
        help="Which space we want to register to. Can be either 'asl' for "
            + "registration to the first volume of the ASL series or "
            + "'structural' for registration to the T1w image. Default "
            + " is 'asl'.",
        default="asl"
    )
    parser.add_argument(
        "--fmap_ap",
        help="Filename for the AP fieldmap for use in distortion correction",
        required=True
    )
    parser.add_argument(
        "--fmap_pa",
        help="Filename for the PA fieldmap for use in distortion correction",
        required=True
    )
    parser.add_argument(
        '--use_t1',
        help="If this flag is provided, the T1 estimates from the satrecov "
            + "will also be registered to ASL-gridded T1 space for use in "
            + "perfusion estimation via oxford_asl.",
        action='store_true'
    )
    parser.add_argument(
        '--sebased',
        help="If this flag is provided, the distortion warps and motion "
            +"estimates will be applied to the MT-corrected but not bias-"
            +"corrected calibration and ASL images. The bias-field will "
            +"then be estimated from the calibration image using HCP's "
            +"SE-based algorithm and applied in subsequent steps.",
        action='store_true'
    )
    parser.add_argument(
        "--mtname",
        help="Filename of the empirically estimated MT-correction"
            + "scaling factors.",
        default=None,
        required="--sebased" in sys.argv
    )
    parser.add_argument(
        "-c",
        "--cores",
        help="Number of cores to use when applying motion correction and "
            +"other potentially multi-core operations. Default is the "
            +f"number of cores your machine has ({mp.cpu_count()}).",
        default=mp.cpu_count(),
        type=int,
        choices=range(1, mp.cpu_count()+1)
    )
    parser.add_argument(
        "--interpolation",
        help="Interpolation order for registrations. This can be any "
            +"integer from 0-5 inclusive. Default is 3. See scipy's "
            +"map_coordinates for more details.",
        default=3,
        type=int,
        choices=range(0, 5+1)
    )
    args = parser.parse_args()
    study_dir = args.study_dir
    sub_id = args.sub_id
    grad_coefficients = args.grads
    target = args.target
    pa_sefm = args.fmap_pa
    ap_sefm = args.fmap_ap
    use_t1 = args.use_t1
    use_sebased = args.sebased
    mt_factors = args.mtname

    # For debug, re-use existing intermediate files 
    force_refresh = True

    # Input, output and intermediate directories
    # Create if they do not already exist. 
    sub_base = op.abspath(op.join(study_dir, sub_id))
    grad_coefficients = op.abspath(grad_coefficients)
    pvs_dir = op.join(sub_base, "T1w", "ASL", "PVEs")
    t1_asl_dir = op.join(sub_base, "T1w", "ASL")
    distcorr_dir = op.join(sub_base, "ASL", "TIs", "DistCorr")
    reg_dir = op.join(sub_base, 'T1w', 'ASL', 'reg')
    t1_dir = op.join(sub_base, "T1w")
    asl_dir = op.join(sub_base, "ASL", "TIs", "STCorr2")
    asl_out_dir = op.join(t1_asl_dir, "TIs", "DistCorr")
    calib_out_dir = op.join(t1_asl_dir, "Calib", "Calib0", "DistCorr") if target=='structural' else op.join(sub_base, "ASL", "Calib", "Calib0", "DistCorr")
    [ os.makedirs(d, exist_ok=True) 
        for d in [pvs_dir, t1_asl_dir, distcorr_dir, reg_dir, 
                  asl_out_dir, calib_out_dir] ]
        
    # Images required for processing 
    asl = op.join(asl_dir, "tis_stcorr.nii.gz")
    struct = op.join(t1_dir, "T1w_acpc_dc_restore.nii.gz")
    struct_brain = op.join(t1_dir, "T1w_acpc_dc_restore_brain.nii.gz")
    struct_brain_mask = op.join(t1_dir, "T1w_acpc_dc_restore_brain_mask.nii.gz")
    asl_vol0 = op.join(asl_dir, "tis_stcorr_vol1.nii.gz")
    if (not op.exists(asl_vol0) or force_refresh) and target=='asl':
        cmd = "fslroi {} {} 0 1".format(asl, asl_vol0)
        sp.run(cmd.split(" "), check=True)

    # Create ASL-gridded version of T1 image 
    t1_asl_grid = op.join(t1_dir, "ASL", "reg", 
                          "ASL_grid_T1w_acpc_dc_restore.nii.gz")
    if (not op.exists(t1_asl_grid) or force_refresh) and target=='asl':
        asl_spc = rt.ImageSpace(asl)
        t1_spc = rt.ImageSpace(struct)
        t1_asl_grid_spc = t1_spc.resize_voxels(asl_spc.vox_size / t1_spc.vox_size)
        nb.save(
            rt.Registration.identity().apply_to_image(struct, t1_asl_grid_spc, order=args.interpolation), 
            t1_asl_grid)
    
    # Create ASL-gridded version of T1 image
    t1_asl_grid_mask = op.join(reg_dir, "ASL_grid_T1w_acpc_dc_restore_brain_mask.nii.gz")
    if (not op.exists(t1_asl_grid_mask) or force_refresh) and target=='asl':
        asl_spc = rt.ImageSpace(asl)
        t1_spc = rt.ImageSpace(struct_brain)
        t1_asl_grid_spc = t1_spc.resize_voxels(asl_spc.vox_size / t1_spc.vox_size)
        t1_mask = binarise_image(struct_brain)
        t1_mask_asl_grid = rt.Registration.identity().apply_to_array(t1_mask, t1_spc, t1_asl_grid_spc, 
                                                                    order=args.interpolation)
        # Re-binarise downsampled mask and save
        t1_asl_grid_mask_array = binary_fill_holes(t1_mask_asl_grid>0.25).astype(np.float32)
        t1_asl_grid_spc.save_image(t1_asl_grid_mask_array, t1_asl_grid_mask) 

    # MCFLIRT ASL using the calibration as reference 
    calib = op.join(sub_base, 'ASL', 'Calib', 'Calib0', 'MTCorr', 'calib0_mtcorr.nii.gz')
    asl = op.join(sub_base, 'ASL', 'TIs', 'tis.nii.gz')
    mcdir = op.join(sub_base, 'ASL', 'TIs', 'MoCo', 'asln2m0.mat')
    asl2calib_mc = rt.MotionCorrection.from_mcflirt(mcdir, asl, calib)

    # Rebase the motion correction to target volume 0 of ASL 
    # The first registration in the series gives us ASL-calibration transform
    calib2asl0 = asl2calib_mc[0].inverse()
    asl_mc = rt.chain(asl2calib_mc, calib2asl0)

    # load the gradient distortion correction warp 
    gdc_path = op.join(sub_base, "ASL", "gradient_unwarp", "fullWarp_abs.nii.gz")
    gdc = rt.NonLinearRegistration.from_fnirt(gdc_path, asl_vol0, 
            asl_vol0, intensity_correct=True, constrain_jac=(0.01,100))

    # get fieldmap names for use with asl_reg
    fmap, fmapmag, fmapmagbrain = [ 
        op.join(sub_base, "ASL", "topup", '{}.nii.gz'.format(s)) 
        for s in [ 'fmap', 'fmapmag', 'fmapmagbrain' ]
    ]

    # get linear registration from asl to structural
    if target == 'asl':
        unreg_img = asl_vol0
    elif target == 'structural':
        # register perfusion-weighted image to structural instead of asl 0
        unreg_img = op.join(sub_base, "ASL", "TIs", "OxfordASL", 
                            "native_space", "perfusion.nii.gz")
    
    # set correct output directory
    distcorr_out_dir = asl_out_dir if target=='structural' else distcorr_dir

    # Initial (linear) asl to structural registration, via first round of asl_reg
    # only need this if target space == structural
    asl2struct_initial_path = op.join(reg_dir, 'asl2struct_final_bbr_fsl.mat')
    if (not op.exists(asl2struct_initial_path) or force_refresh) and target=='structural':
        asl2struct_initial_path_temp = op.join(reg_dir, 'asl2struct_initial_bbr_fsl.mat')
        fsdir = op.join(t1_dir, f'{sub_id}_V1_MR')
        generate_asl2struct_initial(unreg_img, struct, fsdir, reg_dir)
        os.replace(asl2struct_initial_path_temp, asl2struct_initial_path)
    if target == 'structural':
        asl2struct_initial = rt.Registration.from_flirt(asl2struct_initial_path, 
                                                        src=unreg_img, ref=struct)
    elif target == 'asl':
        calib2struct_name = op.join(calib_out_dir, "asl2struct.mat")
        calib2struct = rt.Registration.from_flirt(calib2struct_name,
                                                  calib,
                                                  struct_brain)
        asl2struct_initial = rt.chain(calib2asl0.inverse(), calib2struct)

    # Get brain mask in asl space
    if target == 'asl':
        mask_name = op.join(reg_dir, "asl_vol1_mask_init.nii.gz")
    else:
        mask_name = op.join(reg_dir, "asl_vol1_mask_final.nii.gz")
    if not op.exists(mask_name) or force_refresh:
        asl_mask = generate_asl_mask(struct_brain, unreg_img, asl2struct_initial)
        rt.ImageSpace.save_like(unreg_img, asl_mask, mask_name)

    # get binary WM mask name (using FS' aparc+aseg)
    wmmask = op.join(sub_base, "T1w", "ASL", "reg", "wmmask.nii.gz")

    # load/estimate epi distortion correction
    # if target is asl, we already have a rough estimate from calibration image
    # if target is structural, refine epi distortion correction via asl_reg --finalonly
    if target == 'asl':
        epi_dc_path = op.join(calib_out_dir, "asl2struct_warp.nii.gz")
    else:
        epi_dc_path = op.join(distcorr_dir, "asl2struct_warp_final.nii.gz")
    if (not op.exists(epi_dc_path) and target=='structural') or force_refresh:
        epi_dc_path_temp = op.join(distcorr_dir, 'asl2struct_warp.nii.gz')
        generate_epidc_warp(unreg_img, struct, struct_brain, 
                            mask_name, wmmask, asl2struct_initial, fmap, 
                            fmapmag, fmapmagbrain, distcorr_dir)
        # rename warp so it isn't overwritten
        os.replace(epi_dc_path_temp, epi_dc_path)
    epi_dc = rt.NonLinearRegistration.from_fnirt(epi_dc_path, 
                mask_name, struct, intensity_correct=True, 
                constrain_jac=(0.01,100))

    # if ending in asl space, chain struct2asl transformation
    if target == 'asl':
        struct2calib_reg = op.join(calib_out_dir, "struct2asl.mat")
        struct2calib_reg = rt.Registration.from_flirt(struct2calib_reg,
                                                      src=struct, 
                                                      ref=asl)
        epi_dc = rt.chain(epi_dc, struct2calib_reg)

    # Final ASL transforms: moco, grad dc, 
    # epi dc (incorporating asl->struct reg)
    asl = op.join(sub_base, "ASL", "TIs", "MTCorr", "tis_mtcorr.nii.gz")
    reference = t1_asl_grid if target=='structural' else asl
    if use_sebased and target=='structural':
        asl = op.join(sub_base, "ASL", "TIs", "tis.nii.gz")
        asl = Image(asl)
        asl_sfs = op.join(asl_dir, "combined_scaling_factors.nii.gz")
        asl_sfs = Image(asl_sfs)
        asl = Image(asl.data * asl_sfs.data, header=asl.header)
        mtcorr_name = op.join(sub_base, "ASL", "TIs", "MTCorr", "tis_mtcorr_nobc.nii.gz")
        asl.save(mtcorr_name)
        asl = mtcorr_name
    asl_outpath = op.join(distcorr_out_dir, "tis_distcorr.nii.gz")
    if (not op.exists(asl_outpath) or force_refresh) and target=="structural":
        asl2struct_mc_dc = rt.chain(asl_mc, gdc, epi_dc)
        asl_corrected = asl2struct_mc_dc.apply_to_image(src=asl, 
                                                        ref=reference, 
                                                        cores=args.cores,
                                                        order=args.interpolation)
        nb.save(asl_corrected, asl_outpath)

    # Final calibration transforms: calib->asl, grad dc, 
    # epi dc (incorporating asl->struct reg)
    if use_sebased and target=='structural':
        calib = op.join(sub_base, "ASL", "Calib", "Calib0", "calib0.nii.gz")
        calib = Image(calib)
        mt_sfs = np.loadtxt(mt_factors).reshape(1, 1, -1)
        calib = Image(calib.data * mt_sfs, header=calib.header)
        mtcorr_name = op.join(sub_base, "ASL", "Calib", "Calib0", "MTCorr", "calib0_mtcorr_nobc.nii.gz")
        calib.save(mtcorr_name)
        calib = mtcorr_name
    calib_outpath = op.join(calib_out_dir, "calib0_dcorr.nii.gz")
    if (not op.exists(calib_outpath) or force_refresh) and target=='structural':
        calib2struct_dc = rt.chain(calib2asl0, gdc, epi_dc)
        calib_corrected = calib2struct_dc.apply_to_image(src=calib, 
                                                         ref=reference,
                                                         order=args.interpolation)
        
        nb.save(calib_corrected, calib_outpath)

    # apply distortion corrections to fmapmag.nii.gz
    if use_sebased and target=='structural':
        fmap_reg_dir = op.join(sub_base, "T1w", "ASL", "reg", "fmap")
        bbr_mat = register_fmap(fmapmag, fmapmagbrain, struct, struct_brain, fmap_reg_dir, wmmask)
        fmap2struct_bbr = rt.Registration.from_flirt(bbr_mat, src=fmapmag, ref=struct)
        fmap_struct = fmap2struct_bbr.apply_to_image(src=fmapmag, ref=reference)
        fmap_struct_name = op.join(fmap_reg_dir, "fmapmag_aslstruct.nii.gz")
        nb.save(fmap_struct, fmap_struct_name)

    # Final scaling factors transforms: moco, grad dc, 
    # epi dc (incorporating asl->struct reg)
    sfs_name = op.join(asl_dir, "combined_scaling_factors.nii.gz")
    sfs_outpath = op.join(distcorr_out_dir, "combined_scaling_factors.nii.gz")
    if (not op.exists(sfs_outpath) or force_refresh) and target=="structural":
        # don't chain transformations together if we don't have to
        try:
            asl2struct_mc_dc
        except NameError:
            asl2struct_mc_dc = rt.chain(asl_mc, gdc, epi_dc)
        sfs_corrected = asl2struct_mc_dc.apply_to_image(src=sfs_name, 
                                                        ref=reference, 
                                                        cores=mp.cpu_count())
        nb.save(sfs_corrected, sfs_outpath)
    
    # apply registrations to satrecov-estimated T1 image for use with oxford_asl
    if use_t1:
        est_t1_name = op.join(sub_base, "ASL", "TIs", "SatRecov2", 
                                "spatial", "mean_T1t_filt.nii.gz")
        reg_est_t1_name = op.join(reg_dir, "mean_T1t_filt.nii.gz")
        if (not op.exists(reg_est_t1_name) or force_refresh) and target=='structural':
            asl2struct_dc = rt.chain(asl_mc[0], gdc, epi_dc)
            reg_est_t1 = asl2struct_dc.apply_to_image(src=est_t1_name,
                                                      ref=reference,
                                                      order=args.interpolation)
            nb.save(reg_est_t1, reg_est_t1_name)

    # create ti image in asl space
    slicedt = 0.059
    tis = [1.7, 2.2, 2.7, 3.2, 3.7]
    sliceband = 10
    ti_asl = op.join(sub_base, "ASL", "TIs", "timing_img.nii.gz")
    if (not op.exists(ti_asl) or force_refresh) and target=='asl':
        create_ti_image(asl, tis, sliceband, slicedt, ti_asl)
    
    # transform ti image into t1 space
    ti_t1 = op.join(t1_asl_dir, "timing_img.nii.gz")
    if (not op.exists(ti_t1) or force_refresh) and target=='structural':
        asl2struct = op.join(distcorr_dir, "asl2struct.mat")
        asl2struct = rt.Registration.from_flirt(asl2struct,
                                                src=asl,
                                                ref=struct)
        ti_t1_img = asl2struct.apply_to_image(src=ti_asl,
                                              ref=reference,
                                              order=0,
                                              superfactor=False)
        nb.save(ti_t1_img, ti_t1)

if __name__  == '__main__':

    # study_dir = 'HCP_asl_min_req'
    # sub_number = 'HCA6002236'
    # grad_coefficients = 'HCP_asl_min_req/coeff_AS82_Prisma.grad'
    # sys.argv[1:] = ('%s %s -g %s' % (study_dir, sub_number, grad_coefficients)).split()
    main()
