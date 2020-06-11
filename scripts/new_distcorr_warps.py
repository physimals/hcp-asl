import os
import subprocess as sp
import sys
import os.path as op 
import glob 
import tempfile 

import regtricks as rt
import nibabel as nb
import numpy as np
from fsl.wrappers import fslmaths, bet
from scipy.ndimage import binary_fill_holes

from hcpasl.extract_fs_pvs import extract_fs_pvs
import argparse

def find_field_maps(study_dir, subject_number):
    """
    Find the mbPCASL field maps in the subject's directory.
    The field maps are found in the subject's B session directory. 
    Multiple pairs of field maps are taken in the B session; this 
    function assumes that the mbPCASL field maps are the final 2 
    field map directories in the session.
    """
    scan_dir = op.join(study_dir, subject_number, f'{subject_number}_V1_B', 'scans')
    fm_dirs = sorted(glob.glob(op.join(scan_dir, '*-FieldMap_SE_EPI')))[-2:]
    if op.exists(op.join(fm_dirs[0], 
        f'resources/NIFTI/files/{subject_number}_V1_B_PCASLhr_SpinEchoFieldMap_PA.nii.gz')):
        pa_dir, ap_dir = fm_dirs
    elif op.exists(op.join(fm_dirs[1], 
        f'resources/NIFTI/files/{subject_number}_V1_B_PCASLhr_SpinEchoFieldMap_PA.nii.gz')):
        ap_dir, pa_dir = fm_dirs
    pa_sefm = op.join(pa_dir, f'resources/NIFTI/files/{subject_number}_V1_B_PCASLhr_SpinEchoFieldMap_PA.nii.gz')
    ap_sefm = op.join(ap_dir, f'resources/NIFTI/files/{subject_number}_V1_B_PCASLhr_SpinEchoFieldMap_AP.nii.gz')
    return str(pa_sefm), str(ap_sefm)

def generate_asl2struct_initial(asl, outdir, struct, struct_brain):
    """
    Generate the initial linear transformation between ASL-space and T1w-space
    using asl_reg. This is required as the initalization for the epi distortion 
    correction warp (calculated via asl_reg later on).
    
    Args:
        asl: path to ASL image 
        outdir: path to registration directory, for output 
        struct: path to T1 image, ac_dc_restore
        struct_brain: path to brain-extracted T1 image, ac_dc_restore_brain

    Returns: 
        n/a, file 'asl2struct.mat' will be created in the output dir 
    """
    reg_call = ("asl_reg -i " + asl + " -o " + outdir + " -s " + struct +
                " --sbet=" + struct_brain + " --mainonly")

    sp.run(reg_call.split(), check=True, stderr=sp.PIPE, stdout=sp.PIPE)

def generate_gdc_warp(asl_vol0, coeffs_path, distcorr_dir):
    """
    Generate distortion correction warp via gradient_unwarp. 

    Args: 
        asl_vol0: path to first volume of ASL series
        coeffs_path: path to coefficients file for the scanner (.grad)
        distcorr_dir: directory in which to put output
    
    Returns: 
        n/a, file 'fullWarp_abs.nii.gz' will be created in output dir
    """

    # Need to run in the output directory to make sure files end up in the
    # right place
    pwd = os.getcwd()
    os.chdir(distcorr_dir)
    cmd = ("gradient_unwarp.py {} gdc_corr_vol1.nii.gz siemens -g {}"
            .format(asl_vol0, coeffs_path))
    sp.run(cmd, shell=True)
    os.chdir(pwd)

def generate_wmmask(fsdir, t1_path):
    """ 
    Generate binary WM mask in space of T1 image using FS aparc+aseg

    Args: 
        fsdir: path to freesurfer directory from which aparc+aseg is loaded
        t1_path: path to T1 image (any), for setting output properties

    Returns: 
        np.array logical WM mask in space of T1 image 
    """

    aseg = op.join(fsdir, "mri", "aparc+aseg.mgz")
    aseg_array = nb.load(aseg).get_data()
    wm = np.logical_or(aseg_array == 41, aseg_array == 2)
    wm_t1 = rt.Registration.identity().apply_to_array(wm, aseg, t1_path)
    return wm_t1.round()

def generate_fmaps(pa_ap_sefms, params, config, distcorr_dir): 
    """
    Generate fieldmaps via topup for use with asl_reg. 

    Args: 
        asl_vol0: path to image of stacked blipped images (ie, PEdir as vol0,
            (oPEdir as vol1), in this case stacked as pa then ap)
        params: path to text file for topup --datain, PE directions/times
        config: path to text file for topup --config, other args 
        distcorr_dir: directory in which to put output
    
    Returns: 
        n/a, files 'fmap, fmapmag, fmapmagbrain.nii.gz' will be created in output dir
    """

    # Run topup to get fmap in Hz 
    topup_fmap = op.join(distcorr_dir, 'topup_fmap_hz.nii.gz')        
    cmd = (("topup --imain={} --datain={}".format(pa_ap_sefms, params)
            + " --config={} --out=topup".format(config))
            + " --fout={} --iout={}".format(topup_fmap,
                op.join(distcorr_dir, 'corrected_sefms.nii.gz')))
    sp.run(cmd, shell=True)

    fmap, fmapmag, fmapmagbrain = [ 
        op.join(distcorr_dir, '{}.nii.gz'.format(s)) 
        for s in [ 'fmap', 'fmapmag', 'fmapmagbrain' ]
    ]    

    # Convert fmap from Hz to rad/s
    fmap_spc = rt.ImageSpace(topup_fmap)
    fmap_arr_hz = nb.load(topup_fmap).get_data()
    fmap_arr = fmap_arr_hz * 2 * np.pi
    fmap_spc.save_image(fmap_arr, fmap)

    # Mean across volumes of corrected sefms to get fmapmag
    fmapmag_arr = nb.load(op.join(
                    distcorr_dir, "corrected_sefms.nii.gz")).get_data()
    fmapmag_arr = fmapmag_arr.mean(-1)
    fmap_spc.save_image(fmapmag_arr, fmapmag)

    # Run BET on fmapmag to get brain only version 
    bet(fmap_spc.make_nifti(fmapmag_arr), output=fmapmagbrain)

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

def generate_epidc_warp(asl_vol0_brain, struct, struct_brain, asl_mask,
                       wmmask, asl2struct, fmap, fmapmag, fmapmagbrain, 
                       distcorr_dir):
    """
    Generate EPI distortion correction warp via asl_reg. 

    Args: 
        asl_vol0_brain: path to first volume of ASL series, brain-extracted
        struct: path to T1 image, ac_dc_restore
        struct_brain: path to brain-extracted T1 image, ac_dc_restore_brain
        asl_mask: path to brain mask in ASL space 
        wmmask: path to WM mask in ASL space 
        asl2struct: regtricks.Registration for asl to structural
        fmap: path to topup's field map in rad/s
        fmapmag: path to topup's field map magnitude 
        fmapmagbrain: path to topup's field map magnitude, brain only
        distcorr_dir: path to directory in which to place output 

    Returns: 
        n/a, file 'asl2struct_warp.nii.gz' is created in output directory 
    """

    a2s_fsl = op.join(distcorr_dir, 'asl2struct.mat')
    asl2struct.save_fsl(a2s_fsl, asl_vol0_brain, struct)
    cmd = ("asl_reg -i {} -o {} ".format(asl_vol0_brain, distcorr_dir)
           + "-s {} --sbet={} -m {} ".format(struct, struct_brain, asl_mask)
           + "--tissseg={} --imat={} --finalonly ".format(wmmask, a2s_fsl)
           + "--fmap={} --fmapmag={} ".format(fmap, fmapmag)
           + "--fmapmagbrain={} --pedir=y --echospacing=0.00057 ".format(fmapmagbrain))
    sp.run(cmd, shell=True)




def main(study_dir, sub_id, grad_coefficients):

    # For debug, don't re-use existing intermediate files 
    force_refresh = True

    # Input, output and intermediate directories
    # Create if they do not already exist. 
    sub_base = op.abspath(op.join(study_dir, sub_id))
    pvs_dir = op.join(sub_base, "T1w", "ASL", "PVEs")
    t1_asl_dir = op.join(sub_base, "T1w", "ASL")
    distcorr_dir = op.join(sub_base, "ASL", "TIs", "DistCorr")
    reg_dir = op.join(sub_base, 'T1w', 'ASL', 'reg')
    t1_dir = op.join(sub_base, "T1w")
    asl_dir = op.join(sub_base, "ASL", "TIs", "STCorr", "SecondPass")
    final_out_dir = op.join(t1_asl_dir, "TIs", "DistCorr")
    for d in [pvs_dir, t1_asl_dir, distcorr_dir, reg_dir, final_out_dir]:
        os.makedirs(d, exist_ok=True)

    asl = op.join(asl_dir, "tis_stcorr.nii.gz")
    struct = op.join(t1_dir, "T1w_acpc_dc_restore.nii.gz")
    struct_brain = op.join(t1_dir, "T1w_acpc_dc_restore_brain.nii.gz")
    struct_brain_mask = op.join(t1_dir, "T1w_acpc_dc_restore_brain_mask.nii.gz")
    asl_vol0 = op.join(asl_dir, "tis_stcorr_vol1.nii.gz")
    if not op.exists(asl_vol0) or force_refresh:
        cmd = "fslroi {} {} 0 1".format(asl, asl_vol0)

    # Create ASL-gridded version of T1 image 
    t1_asl_grid = op.join(t1_dir, "ASL", "reg", "ASL_grid_T1w_acpc_dc_restore.nii.gz")
    if not op.exists(t1_asl_grid) or force_refresh:
        asl_spc = rt.ImageSpace(asl)
        t1_spc = rt.ImageSpace(struct)
        t1_asl_grid_spc = t1_spc.resize_voxels(asl_spc.vox_size / t1_spc.vox_size)
        nb.save(rt.Registration.identity().apply_to_image(struct, t1_asl_grid_spc), t1_asl_grid)

    # MCFLIRT ASL using the calibration as reference 
    calib = op.join(sub_base, 'ASL', 'Calib', 'Calib0', 'calib0.nii.gz')
    asl = op.join(sub_base, 'ASL', 'TIs', 'tis.nii.gz')
    mcdir = op.join(sub_base, 'ASL', 'TIs', 'MoCo', 'asl2calib.mat')
    if not op.exists(mcdir) or force_refresh:
        asl2calib_mc = rt.mcflirt(asl, reffile=calib)
        asl2calib_mc.save_fsl(mcdir, asl, calib)
    else: 
        asl2calib_mc = rt.MotionCorrection.from_mcflirt(mcdir, asl, calib)

    # Rebase the motion correction to target volume 0 of ASL 
    # The first registration in the series gives us ASL-calibration transform
    calib2asl0 = asl2calib_mc[0].inverse()
    asl_mc = rt.chain(asl2calib_mc, calib2asl0)

    # Generate the gradient distortion correction warp 
    gdc_path = op.join(distcorr_dir, 'fullWarp_abs.nii.gz')
    if not op.exists(gdc_path) or force_refresh:
        gdc = generate_gdc_warp(asl_vol0, grad_coefficients, distcorr_dir)
    gdc = rt.NonLinearRegistration.from_fnirt(gdc_path, asl_vol0, 
            asl_vol0, intensity_correct=True, constrain_jac=True)

    # Stack the cblipped images together for use with topup 
    pa_sefm, ap_sefm = find_field_maps(study_dir, sub_number)
    pa_ap_sefms = op.join(distcorr_dir, 'merged_sefms.nii.gz')
    if not op.exists(pa_ap_sefms) or force_refresh:
        rt.ImageSpace.save_like(pa_sefm, np.stack((
                                nb.load(pa_sefm).get_data(), 
                                nb.load(ap_sefm).get_data()), 
                                axis=-1), 
                                pa_ap_sefms)
    topup_params = op.join(distcorr_dir, 'topup_params.txt')
    topup_config = "b02b0.cnf"  # Note this file doesn't exist in scope, 
                                # but topup knows where to find it 

    # Generate fieldmaps for use with asl_reg (via topup)
    fmap, fmapmag, fmapmagbrain = [ 
        op.join(distcorr_dir, '{}.nii.gz'.format(s)) 
        for s in [ 'fmap', 'fmapmag', 'fmapmagbrain' ]
    ]          
    if not all([ op.exists(p) for p in [fmap, fmapmag, fmapmagbrain] ]):# or force_refresh:                  
        generate_fmaps(pa_ap_sefms, topup_params, topup_config, distcorr_dir)

    # Initial (linear) asl to structural registration, via first round of asl_reg
    asl2struct_initial_path = op.join(reg_dir, 'asl2struct.mat')
    if not op.exists(asl2struct_initial_path) or force_refresh:
        generate_asl2struct_initial(asl_vol0, reg_dir, struct, struct_brain)
    asl2struct_initial = rt.Registration.from_flirt(asl2struct_initial_path, asl, struct)

    # Get brain mask in asl space
    asl_mask_path = op.join(reg_dir, "asl_vol1_mask.nii.gz")
    if not op.exists(asl_mask_path) or force_refresh:
        asl_mask = generate_asl_mask(struct_brain, asl_vol0, asl2struct_initial)
        rt.ImageSpace.save_like(asl, asl_mask, asl_mask_path)
        asl_mask = asl_mask_path

    # Brain extract volume 0 of asl series 
    asl_vol0_brain = op.join(sub_base, "ASL", "TIs", "STCorr", "SecondPass", "tis_stcorr_vol1.nii.gz")
    if not op.exists(asl_vol0_brain) or force_refresh:
        bet(asl_vol0, asl_vol0_brain)

    # # Generate PVEs in T1 ASL-gridded space 
    # pve_base = op.join(sub_base, "T1w", "ASL", "PVEs", "pve")
    # if not glob.glob(pve_base + '*') or force_refresh: 
    #     aparc_aseg = op.join(sub_base, "T1w", "aparc+aseg.nii.gz")
    #     gen_pves(op.dirname(struct), t1_asl_grid, pve_base)

    # Generate a binary WM mask in the space of the T1 (using FS' aparc+aseg)
    wmmask = op.join(sub_base, "T1w", "wmmask.nii.gz")
    if not op.exists(wmmask) or force_refresh:
        fsdir = op.join(t1_dir, "%s_V1_MR" % sub_id)
        wmmask_img = generate_wmmask(fsdir, struct)
        rt.ImageSpace.save_like(struct, wmmask_img, wmmask)

    # Generate the EPI distortion correction warp via asl_reg --final
    epi_dc_path = op.join(distcorr_dir, 'asl2struct_warp.nii.gz')
    if not op.exists(epi_dc_path):# or force_refresh:
        generate_epidc_warp(asl_vol0_brain, struct, struct_brain, asl_mask_path,
                                wmmask, asl2struct_initial, fmap, fmapmag, fmapmagbrain, 
                                distcorr_dir)
    epi_dc = rt.NonLinearRegistration.from_fnirt(epi_dc_path, 
                asl_mask_path, struct, intensity_correct=True, constrain_jac=True)

    # Stick them all together!
    asl2asl_mc_dc_chain = rt.chain(asl_mc, gdc, epi_dc)
    asl_corrected = asl2asl_mc_dc_chain.apply_to_image(src=asl, ref=t1_asl_grid, cores=8)
    outpath = op.join(final_out_dir, "tis_distcorr_rt.nii.gz")
    nb.save(asl_corrected, outpath)

if __name__  == '__main__':

    study_dir = 'HCP_asl_min_req'
    sub_number = 'HCA6002236'
    grad_coefficients = op.abspath('HCP_asl_min_req/coeff_AS82_Prisma.grad')

    sys.argv[1:] = ('%s %s -g %s' % (study_dir, sub_number, grad_coefficients)).split()

    # argument handling
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "study_dir",
        help="Path of the base study directory."
    )
    parser.add_argument(
        "sub_number",
        help="Subject number."
    )
    parser.add_argument(
        "-g",
        "--grads",
        help="Filename of the gradient coefficients for gradient"
            + "distortion correction (optional)."
    )

    args = parser.parse_args()
    main(args.study_dir, args.sub_number, args.grads)
