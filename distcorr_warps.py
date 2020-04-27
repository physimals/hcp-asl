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
    
def calc_fmaps(pa_ap_sefms, pars_filepath, cnf_file, distcorr_dir, out_basename, topup_fmap, 
                fmap_rads, fmapmag, fmapmagbrain):

    """
    Use topup to generate a fieldmap for distortion correction in asl_reg 
    after conversion from Hz to rads/s. Topup also outputs the SEFMs 
    corrected for EPI distortion, which are then available for use as the 
    fieldmap magnitude image and  brain-extracted fieldmap mag image.
    """
    
    topup_call = ("topup --imain=" + distcorr_dir + "/" + pa_ap_sefms + 
                    " --datain=" + pars_filepath + " --config=" + cnf_file + 
                    " --out=" + distcorr_dir + "/" + out_basename + " --fout=" +
                    distcorr_dir + "/" + out_fmap + " --iout=" + distcorr_dir +
                    "corrected_sefms.nii.gz")
    
    print(topup_call)
    # sp.run(topup_call.split(), check=True, stderr=PIPE, stdout=sp.PIPE)

    convert_torads = ("fslmaths " + distcorr_dir + "/" + out_fmap + 
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
                        moco_xfms, calib_orig, calib_T1space, calib_xfms):
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
    print(asl_apply_call)
    print(asl_jaco_call)
    print(calib_apply_call)
    print(calib_jaco_call)

    # sp.run(asl_apply_call.split(), check=True, stderr=PIPE, stdout=PIPE)
    # sp.run(asl_jaco_call.split(), check=True, stderr=PIPE, stdout=PIPE)
    # sp.run(calib_apply_call.split(), check=True, stderr=PIPE, stdout=PIPE)
    # sp.run(calib_jaco_call.split(), check=True, stderr=PIPE, stdout=PIPE)

if __name__ == "__main__":
    # fill in function calls

    temp_dir = something
    changedir = ("cd " + temp_dir)
    sp.run(changedir.split())

    # Call function to generate GDC warp:
    asl1 = 
    grad_coeffs = 
    oph 

