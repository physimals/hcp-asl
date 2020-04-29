#! /usr/bin/env python3
""" 
 Script for gradient and EPI ddistortion correction of HCP ASL data

 Script to generate gradient distortion correction, and EPI distortion
 correction warps for the HCP ASL data, then apply the warps along 
 with motion correction to the ASL-gridded T1w-space, so that all 
 transformations are done in a single interpolation step.

 F. A. Kennedy McConnell - April 2020
"""

import nibabel as nb
import numpy as np
import os
import subprocess as sp
import regtools as rt
import tempfile as td
import shutil

from extract_fs_pvs import extract_fs_pvs

# Generate gradient distortion correction warp
def calc_gdc_warp(asldata_vol1, coeffs_loc, oph):
    """
    Generate warp for gradient distortion correction using siemens 
    coefficients file.
    """

    gdc_call = ("python gradient_unwarp.py " + asldata_vol1 + " gdc_corr_vol1.nii.gz " +
                    "siemens -g " + coeffs_loc)

    print(gdc_call)
    #sp.run(gdc_call.split(), check=True, stderr=sp.PIPE, stdout=sp.PIPE)

    gdc_warp_call = ("convertwarp --abs --ref=" + oph + "/gdc_corr_vol1.nii.gz " +
                    "--warp1" + oph + "/fullWarp_abs.nii.gz --relout --out=" + 
                    oph + "/gdc_warp.nii.gz")

    print(gdc_warp_call)
    #sp.run(gdc_warp_call.split(), check=True, stderr=sp.PIPE, stdout=sp.PIPE)

def produce_topup_params(pars_filepath):
    """
    Generate a file containing the parameters used by topup to generate EPI distortion
    correction fieldmap and output this into subject's directory. -- This could be a 
    global/config file if that is preferable
    """
    with open(pars_filepath, "a") as t_pars:
        t_pars.write("0 1 0 0.04845" + "\n")
        t_pars.write("0 -1 0 0.04845" + "\n")
    
def calc_fmaps(pa_sefm, ap_sefm, pa_ap_sefms, pars_filepath, cnf_file, distcorr_dir, out_basename, topup_fmap, 
                fmap_rads, fmapmag, fmapmagbrain):

    """
    Use topup to generate a fieldmap for distortion correction in asl_reg 
    after conversion from Hz to rads/s. Topup also outputs the SEFMs 
    corrected for EPI distortion, which are then available for use as the 
    fieldmap magnitude image and  brain-extracted fieldmap mag image.
    """
    merge_sefm_call = ("fslmerge -t " + pa_ap_sefms + " " + pa_sefm + " " + ap_sefm)
    print(merge_sefm_call)
    # sp.run(merge_sefm_call.split(), check=True, stderr=PIPE, stdout=PIPE)
    
    topup_call = ("topup --imain=" + distcorr_dir + "/" + pa_ap_sefms + 
                    " --datain=" + pars_filepath + " --config=" + cnf_file + 
                    " --out=" + distcorr_dir + "/" + out_basename + " --fout=" +
                    distcorr_dir + "/" + topup_fmap + " --iout=" + distcorr_dir +
                    "corrected_sefms.nii.gz")
    
    print(topup_call)
    # sp.run(topup_call.split(), check=True, stderr=PIPE, stdout=sp.PIPE)

    convert_torads = ("fslmaths " + distcorr_dir + "/" + topup_fmap + 
                        " -mul 3.14159 -mul 2 " + distcorr_dir + "/" + fmap_rads)
    print(convert_torads)
    # sp.run(convert_torads.split(), check=True, stderr=PIPE, stdout=PIPE)

    mean_fmapmag_call = ("fslmaths " + distcorr_dir + "/corrected_sefms.nii.gz -Tmean " + 
                            distcorr_dir + "/" + fmapmag)
    print(mean_fmapmag_call)
    # sp.run(mean_fmapmag_call.split(), check=True, stderr=PIPE, stdout=PIPE)

    bet_fmapmag_call = ("bet " + distcorr_dir + "/" + fmapmag + " " + distcorr_dir +
                            "/" + fmapmagbrain)
    print(bet_fmapmag_call)
    # sp.run(bet_fmapmag_call.split(), check=True, stderr=PIPE, stdout=PIPE)

def gen_initial_trans(regfrom, outdir, struct, struct_brain):
    """
    Generate the initial linear transformation between ASL-space and T1w-space
    using asl_reg. This is required as the initalization for the registration
    which generates the distortion correcting warp.
    """
    reg_call = ("asl_reg -i " + regfrom + " -o " + outdir + " -s " + struct +
                " --sbet=" + struct_brain + " --mainonly")
    print(reg_call)
    # sp.run(reg_call.split(), check=True, stderr=PIPE, stdout=PIPE)

def gen_asl_mask(struct_brain, struct_bet_mask, regfrom, asl2struct, asl_mask,
                struct2asl):
    """
    Generate a mask of the brain in the space of the first ASL volume. The
    make is required for use in asl_reg when it is called for the purpose
    of generating te distortion correction warp.
    """
    invert_reg = ("convert_xfm -omat " + struct2asl + " -inverse " + asl2struct)

    sbrain_call = ("fslmaths " + struct_brain + " -bin " + struct_bet_mask)
    trans_call = ("flirt -in " + struct_bet_mask + " -ref " + regfrom + " -applyxfm -init " + 
                    struct2asl + " -out " + asl_mask + " -interp trilinear -paddingsize 1")
    fill_call = ("fslmaths " + asl_mask + " -thr 0.25 -bin -fillh " + asl_mask)
    hdr_call = ("fslcpgeom " + regfrom + " " + asl_mask)

    print(invert_reg)
    print(sbrain_call)
    print(trans_call)
    print(fill_call)
    print(hdr_call)

    # sp.run(invert_reg.split(), check=True, stderr=PIPE, stdout=PIPE)
    # sp.run(sbrain_call.split(), check=True, stderr=PIPE, stdout=PIPE)
    # sp.run(trans_call.split(), check=True, stderr=PIPE, stdout=PIPE)
    # sp.run(fill_call.split(), check=True, stderr=PIPE, stdout=PIPE)
    # sp.run(hdr_call.split(), check=True, stderr=PIPE, stdout=PIPE)

def gen_pves(aparc_aseg, t1, asl, fileroot):
    """
    Generate partial volume estimates from freesurfer segmentations of the cortex
    and subcortical structures.
    """    

    # FIXME: superfactor and cores are at debug settings here. 
    pvs_stacked = extract_fs_pvs(aparc_aseg, t1, asl, superfactor=2, cores=1)
    hdr = pvs_stacked.header 
    aff = pvs_stacked.affine 
    for idx, suffix in enumerate(['GM', 'WM', 'CSF']):
        p = "{}_{}".format(fileroot, suffix)
        nii = nb.Nifti2Image(pvs_stacked.dataobj[...,idx], aff, header=hdr)
        nb.save(nii, p)
    

def calc_distcorr_warp(regfrom, distcorr_dir, struct, struct_brain, mask, tissseg,
                        asl2struct_trans, fmap_rads, fmapmag, fmapmagbrain, asl_grid_T1,
                        gdc_warp):
    """
    Warp is generated from asl_reg, which performs EPI distortion correction of 
    ASL data while transforming the data to alignment with T1w-space. This warp is
    then merged with the gradient distortion correction warp if that is available.
    """

    ## Use asl_reg to inital EPI distcorr warp
    reg_call = ("asl_reg -i " + regfrom + " -o " + distcorr_dir + " -s " +
                struct + " --sbet=" + struct_brain + " -m " + mask + " --tissseg " +
                tissseg + " --imat " + asl2struct_trans + " --finalonly --fmap=" + 
                fmap_rads + " --fmapmag=" + fmapmag + " --fmapmagbrain=" + fmapmagbrain +
                " --pedir=y --echospacing=0.00057")
    print(reg_call)
    # sp.run(reg_call.split(), check=True, stderr=PIPE, stdout=PIPE)
    
    # Add the gradient distortion correction warp to the EPI distortion correction warp
    if os.path.isfile(gdc_warp):
        merge_warp_call = ("convertwarp -r " + asl_grid_T1 + " -o " + distcorr_dir + 
                        "/distcorr_warp -w " + distcorr_dir + "/asl2struct_warp --warp2=" + 
                        gdc_warp + " --rel")
        print(merge_warp_call)
        # sp.run(merge_warp_call.split(), check=True, stderr=PIPE, stdout=PIPE)
    else:
        print("Gradtient distortion correction not applied")

# calculate the jacobian of the warp for intensity correction
def calc_warp_jacobian(distcorr_dir):
    """
    Calculation of the Jacobian of the combined distortion correction for subsequent 
    Jacobian intensity scaling
    """
    utils_call1 = ("fnirtfileutils -i " + distcorr_dir + "/distcorr_warp -f spline -o " + 
                    distcorr_dir + "/distcorr_warp_coeff")
    utils_call2 = ("fnirtfileutils -i " + distcorr_dir + "/distcorr_warp_coeff -j " +
                    distcorr_dir + "/distcorr_jacobian")
    hdr_call = ("fslcpgeom " + distcorr_dir + "/distcorr_warp " + distcorr_dir + 
                    "/distcorr_jacobian")

    print(utils_call1)
    print(utils_call2)
    print(hdr_call)

    # sp.run(utils_call1.split(), check=True, stderr=PIPE, stdout=PIPE)
    # sp.run(utils_call2.split(), check=True, stderr=PIPE, stdout=PIPE)
    # sp.run(hdr_call.split(), check=True, stderr=PIPE, stdout=PIPE)

# apply the combined distortion correction warp
def apply_distcorr_warp(asldata_orig, T1space_ref, asldata_T1space, distcorr_dir,
                        moco_xfms, calib_orig, calib_T1space, calib_xfms, sfacs_orig,
                        sfacs_T1space):
    """
    Application of motion correction, and EPI and gradient distortion corrections, 
    and transformation to T1w-space in one step. Leaving the ASL and calibration data 
    in ASL-gridded T1w-space.
    """

    asl_apply_call = ("applywarp -i " + asldata_orig + " -r " + T1space_ref + " -o " +
                    asldata_T1space + " --premat=" + moco_xfms + " -w " + 
                    distcorr_dir + "/distcorr_warp" + " --rel --interp=trilinear" +
                    " --paddingsize=1 --super --superlevel=a")

    asl_jaco_call = ("fslmaths " + asldata_T1space + " -mul " + distcorr_dir + 
                    "/distcorr_jacobian " + asldata_T1space)

    calib_apply_call = ("applywarp -i " + calib_orig + " -r " + T1space_ref + " -o " +
                    calib_T1space + " --premat=" + calib_xfms + " -w " + 
                    distcorr_dir + "/distcorr_warp" + " --rel --interp=trilinear" +
                    " --paddingsize=1 --super --superlevel=a")

    calib_jaco_call = ("fslmaths " + calib_T1space + " -mul " + distcorr_dir + 
                    "/distcorr_jacobian " + calib_T1space)

    sfacs_apply_call = ("applywarp -i " + sfacs_orig + " -r " + T1space_ref + " -o " +
                    sfacs_T1space + " --premat=" + moco_xfms + " -w " + 
                    distcorr_dir + "/distcorr_warp" + " --rel --interp=trilinear" +
                    " --paddingsize=1 --super --superlevel=a")

    sfacs_jaco_call = ("fslmaths " + sfacs_T1space + " -mul " + distcorr_dir + 
                    "/distcorr_jacobian " + sfacs_T1space)
    print(asl_apply_call)
    print(asl_jaco_call)
    print(calib_apply_call)
    print(calib_jaco_call)
    print(sfacs_apply_call)
    print(sfacs_jaco_call)

    # sp.run(asl_apply_call.split(), check=True, stderr=PIPE, stdout=PIPE)
    # sp.run(asl_jaco_call.split(), check=True, stderr=PIPE, stdout=PIPE)
    # sp.run(calib_apply_call.split(), check=True, stderr=PIPE, stdout=PIPE)
    # sp.run(calib_jaco_call.split(), check=True, stderr=PIPE, stdout=PIPE)
    # sp.run(sfacs_apply_call.split(), check=True, stderr=PIPE, stdout=PIPE)
    # sp.run(sfacs_jaco_call.split(), check=True, stderr=PIPE, stdout=PIPE)

if __name__ == "__main__":
    # fill in function calls

    oph = (study_dir + "/" + sub_num + "ASL/TIs/DistCorr")
    changedir = ("cd " + oph)
    sp.run(changedir.split())

    # Generate ASL-gridded T1-aligned T1w image for use as a reg reference
    t1 = (study_dir + "/" + sub_num + "/T1w/T1w_acpc_dc_restore.nii.gz")
    t1_brain = (study_dir + "/" + sub_num + "/T1w/T1w_acpc_dc_restore_brain.nii.gz")

    asl = (study_dir + "/" + sub_num + "/ASL/TIs/STCorr/Second_pass/tis_stcorr.nii.gz") 
    t1_asl_res = (study_dir + "/" + sub_num + "/T1w/ASL/reg/ASL_grid_T1w_acpc_dc_restore.nii.gz")

    asl_v1 = (study_dir + "/" + sub_num + "/ASL/TIs/STCorr/Second_pass/tis_stcorr_vol1.nii.gz")
    first_asl_call = ("fslroi " + asl + " " + asl_v1 + " 0 1")
    print(first_asl_call)
    # sp.run(first_asl_call.split(), check=True, stderr=PIPE, stdout=PIPE)

    t1_spc = rt.ImageSpace(t1)
    asl_spc = rt.ImageSpace(asl_v1)
    t1_spc_asl = t1_spc.resize_voxels(asl_spc.voxel_size / t1_spc.voxel_size)
    r = rt.Registration.identity()
    t1_asl = r.apply_to_image(t1, t1_spc_asl)
    nb.save(t1_asl, t1_asl_res)

    # Check .grad coefficients are available and call function to generate 
    # GDC warp if they are:
    #### Note - I've stored the .grad file in "{$HCPPIPEDIR_Config}", which is the directory ####
    #### listed below. This might not be the right place for this file.                      ####
    grad_coeffs = os.path.expanduser("~/projects/Pipelines/global/config/coeff_AS82_Prisma.grad")
    if os.path.isfile(grad_coeffs):
        calc_gdc_warp(asl_v1, grad_coeffs, oph)
    else:
        print("Gradient coefficients not available")
     
    # output file of topup parameters to subject's distortion correction dir
    pars_filepath = (oph + "/topup_params.txt")
    produce_topup_params(pars_filepath)

    # generate EPI distortion correction fieldmaps for use in asl_reg
    pa_sefm = (sub_num + "/" + sub_num + "_V1_A/30_Fieldmap_SE_EPI/NIFTI/HCA6002236_V1_B_PCASLhr_SpinEchoFieldMap_PA.nii.gz")
       #/mnt/hgfs/Postdoc_data_files/HCP/Aging/HCA6002236/B_scans/30_FieldMap_SE_EPI/NIFTI
    ap_sefm = (sub_num + "/" + sub_num + "_V1_A/30_Fieldmap_SE_EPI/NIFTI/HCA6002236_V1_B_PCASLhr_SpinEchoFieldMap_AP.nii.gz")
    pa_ap_sefms = (oph + "/merged_sefms.nii.gz")
    cnf_file = "b02b0.cnf"
    out_basename = (oph + "/topup_result")
    topup_fmap = (oph + "/topup_result_fmap_hz.nii.gz")
    fmap_rads = (oph + "/fmap_rads.nii.gz")
    fmapmag = (oph + "/fmapmag.nii.gz")
    fmapmagbrain = (oph + "/fmapmag_brain.nii.gz")
    calc_fmaps(pa_sefm, ap_sefm, pa_ap_sefms, pars_filepath, cnf_file, oph, out_basename, topup_fmap, 
                fmap_rads, fmapmag, fmapmagbrain)
    
    # Calculate initial linear transformation from ASL-space to T1w-space
    asl_v1_brain = (study_dir + "/" + sub_num + "/ASL/TIs/STCorr/Second_pass/tis_stcorr_vol1_brain.nii.gz")
    bet_regfrom_call = ("bet " + asl_v1 + " " + asl_v1_brain)
    print(bet_regfrom_call)
    # sp.run(bet_regfrom_call.split(), check=True, stderr=PIPE, stdout=PIPE)

    outdir = (study_dir + "/" + sub_num + "/T1w/ASL/reg")
    gen_initial_trans(asl_v1_brain, outdir, t1, t1_brain)

    # Generate a brain mask in the space of the 1st ASL volume
    asl2struct = (outdir + "/asl2struct.mat")
    t1_brain_mask = (outdir + "/T1w_acpc_dc_restore_brain_mask.nii.gz")
    asl_mask = (outdir + "/asl_vol1_mask.nii.gz")
    struct2asl = (outdir + "/struct2asl.mat")

    gen_asl_mask(t1_brain, t1_brain_mask, asl_v1_brain, asl2struct, asl_mask,
                struct2asl)

    # Generate PVEs
    aparc_aseg = (study_dir + "/" + sub_num + "/T1w/aparc+aseg.nii.gz")
    pve_files = (study_dir + "/" + sub_num + "/T1w/ASL/PVEs/pve")
    gen_pves(aparc_aseg, t1, asl, pve_files)
    

    # Calculate the overall distortion correction warp
    tissseg = 
    asl2str_trans = (outdir + "/asl2struct.mat")
    gdc_warp = (oph + "/gdc_warp.nii.gz")
    calc_distcorr_warp(asl_v1_brain, oph, t1, t1_brain, asl_mask, tissseg,
                        asl2str_trans, fmap_rads, fmapmag, fmapmagbrain, t1_asl_res,
                        gdc_warp)

    # Calculate the Jacobian of the distortion correction warp 
    calc_warp_jacobian(oph)

    # apply the combined distortion correction warp with motion correction
    # to move asl data, calibrationn images, and scaling factors into 
    # ASL-gridded T1w-aligned space
    
    T1w_oph = (study_dir + "/" + sub_num + "/T1w/ASL/TIs/DistCorr")

    asl_distcorr = (T1w_oph + "/tis_distcorr.nii.gz")
    moco_xfms = (study_dir + "/" + sub_num + "/ASL/TIs/MoCo/asln2asl0.mat") #will this work?

    # only correcting and transforming the 1st of the calibration images at the moment
    calib_orig = (study_dir + "/" + sub_num + "/ASL/Calib/Calib0/MTCorr/calib0_mtcorr.nii.gz")
    calib_distcorr = (study_dir + "/" + sub_num + "/T1w/ASL/Calib/Calib0/DistCorr/calib0_dcorr.nii.gz")
    calib_inv_xfm = (study_dir + "/" + sub_num + "/ASL/TIs/MoCo/asln2m0.mat/MAT_0000")
    calib_xfm = (study_dir + "/" + sub_num + "/ASL/TIs/MoCo/calibTOasl1.mat")

    sfacs_orig = (study_dir + "/" + sub_num + "/ASL/TIs/STCorr/Second_pass/st_scaling_factors.nii.gz")
    sfacs_distcorr = (T1w_oph + "/st_scaling_factors.nii.gz")

    invert_call = ("convert_xfm -omat " + calib_xfm + " -inverse " + calib_inv_xfm)
    print(invert_call)
    # sp.run(invert_call.split(), check=True, stderr=PIPE, stdout=PIPE)

    apply_distcorr_warp(asl, t1_asl_res, asl_distcorr, T1w_oph,
                        moco_xfms, calib_orig, calib_distcorr, calib_xfm, sfacs_orig,
                        sfacs_distcorr)


