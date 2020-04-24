# Script for gradient and EPI ddistortion correction of HCP ASL data

### Script to generate gradient distortion correction, and EPI distortion
### correction warps for the HCP ASL data, then apply the warps along 
### with motion correction to the ASL-gridded T1w-space, so that all 
### transformations are done in a single interpolation step.

import nibabel as nb
import numpy as np
import os
import subprocess as sp

# Generate gradient distortion correction warp
def calc_gdc_warp(asldata_vol1, coeffs_loc):
    ##### TO DO #####
    # Got to be careful with the command below - gradient_unwarp.py 
    # will output files directly into the current working directory...
    gdc_call = ("python gradient_unwarp.py " + asldata_vol1 + " gdc_corr_vol1.nii.gz " +
                    "siemens -g " + coeffs_loc)

    print(gdc_call)
    #sp.run(gdc_call.split(), check=True, stderr=sp.PIPE, stdout=sp.PIPE)

    gdc_warp_call = ("convertwarp --abs --ref=gdc_corr_vol1.nii.gz " +
                    "--warp1=fullWarp_abs.nii.gz --relout --out=gdc_warp_path")

    print(gdc_warp_call)
    #sp.run(gdc_warp_call.split(), check=True, stderr=sp.PIPE, stdout=sp.PIPE)

def produce_topup_params(pars_filepath):
    with open(pars_filepath, "a") as t_pars:
        t_pars.write("0 1 0 0.04845" + "\n")
        t_pars.write("0 -1 0 0.04845" + "\n")
    
def calc_fmaps(pa_ap_sefms, pars_filepath, cnf_file, out_basename, out_fmap, out_corrected_sefms, fmapmag, fmapmagbrain):
    
    topup_call = ("topup --imain=" + pa_ap_sefms + " --datain=" + pars_filepath +
                    " --config=" + cnf_file + " --out=" + out_basename +
                    " --fout=" + out_fmap + " --iout=out_corrected_sefms")
    
    print(topup_call)
    # sp.run(topup_call.split(), check=True, stderr=PIPE, stdout=sp.PIPE)

    convert_torads = ("fslmaths " + out_fmap + " -mul 3.14159 -mul 2 " 
                        + fmap_rads)
    print(convert_torads)
    # sp.run(convert_torads.split(), check=True, stderr=PIPE, stdout=PIPE)

    mean_fmapmag_call = ("fslmaths " + out_corrected_sefms + " -Tmean " + 
                            fmapmag)
    print(mean_fmapmag_call)
    # sp.run(mean_fmapmag_call.split(), check=True, stderr=PIPE, stdout=PIPE)

    bet_fmapmag_call = ("bet " + fmapmag + " " + fmapmagbrain)
    print(bet_fmapmag_call)
    # sp.run(bet_fmapmag_call.split(), check=True, stderr=PIPE, stdout=PIPE)

def calc_distcorr_warp(regfrom, distcorr_dir, ):
    print("Performing distortion correction using asl_reg")

    ## Use asl_reg to inital EPI distcorr warp
    reg_call = ("asl_reg -i " + regfrom + " -o " + distcorr_dir + " -s " +
                struct + " --sbet=" + struct_brain + " -m " + mask + " --tissseg " +
                tissseg + " --imat " + asl2struct_trans " --finalonly --fmap=" + 
                fmap_rads + " --fmapmag=" + fmapmag + " --fmapmagbrain=" + fmapmagbrain +
                " --pedir=y --echospacing=0.00057")
    print(reg_call)
    # sp.run(reg_call.split(), check=True, stderr=PIPE, stdout=PIPE)
    
    # get the warp into ASL space
    convertwarp_call = ("convertwarp -r " + regfrom + " -o " + distcorr_dir + "/epi_dist_warp" + 
                        " -w " + distcorr_dir +   "/asl2struct_warp.nii.gz --postmat=" + 
                        struct2asl_trans + " --rel")
    print(convertwarp_call)
    # sp.run(convertwarp_call.split(), check=True, stderr=PIPE, stdout=PIPE)

    # Add the gradient distortion correction warp to the EPI distortion correction warp
    ##### Need to make this conditional on the existence of the GDC warp
    merge_warp_call = ("convertwarp -r " + regfrom + " -o " + distcorr_dir + 
                        "/distcorr_warp -w " + distcorr_dir + "/epi_dist_warp --warp2=" + 
                        gdc_warp + " --rel")
    print(merge_warp_call)
    # sp.run(merge_warp_call.split(), check=True, stderr=PIPE, stdout=PIPE)

# calculate the jacobian of the warp for intensity correction
def calc_warp_jacobian():
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
def apply_distcorr_warp(distcorr_dir, motion_transforms):

    apply_call = ("applywarp -i " + asldata_orig " -r " + T1space_ref + " -o " +
                    asldata_T1space + " --premat=" + moco_xfms + " -w " + 
                    distcorr_dir + "/distcorr_warp" + " --rel --interp=trilinear" +
                    " --paddingsize=1 --super --superlevel=a")

    ### Not sure how to do the ASL data Jacobian intensity scaling given that 
    ### the data is being transformed into ASL-gridded T1w space - what are the 
    ### properties of the Jacobian? The warp will be resolution independent, but
    ### obviously isn't space-independent.

    # If I were trying to stay in ASL space then the bash command could just be:
    # fslmaths $tempdir/asldata -mul $jacobian $tempdir/asldata
    # then it would redo the averaging of the asl data to reproduce regfrom - 
    # which wouldnt't be needed here...

if __name__ == "__main__":
    # fill in function calls