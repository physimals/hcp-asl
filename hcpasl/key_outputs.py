from shutil import copy
import subprocess


def copy_key_outputs(path):

    source_path_T1 = path + "/T1w/ASL/TIs/OxfordASL/native_space/"
    destination_path_T1 = path + "/T1w/ASL/"

    source_path_MNI = path + "/MNI/ASL/Results/OutputtoCIFTI/"
    source_path_MNI_pv = path + "/MNI/ASL/Results/pvcorr/OutputtoCIFTI/"
    destination_path_MNI = path + "/MNINonLinear/ASL/"
    pv_prefix = "pvcorr"

    nonpv_img_files = ["perfusion_calib.nii.gz", \
               "perfusion_var_calib.nii.gz", \
               "arrival.nii.gz", \
               "arrival_var.nii.gz", \
               "aCBV_calib.nii.gz"]

    nonpv_txt_files = ["perfusion_calib_gm_mean.txt", \
               "perfusion_calib_wm_mean.txt", \
               "arrival_gm_mean.txt", \
               "arrival_wm_mean.txt"]

    pv_img_files = ["perfusion_calib_masked.nii.gz", \
            "perfusion_val_calib_masked.nii.gz", \
            "perfusion_wm_calib_masked.nii.gz", \
            "perfusion_wm_var_calib_masked.nii.gz", \
            "arrival_masked.nii.gz", \
            "arrival_var_masked.nii.gz", \
            "arrival_wm_masked.nii.gz", \
            "arrival_wm_var_masked.nii.gz", \
            "aCBV_calib.nii.gz"]

    pv_txt_files = ["perfusion_calib_gm_mean.txt", \
                "perfusion_wm_calib_wm_mean.txt", \
                "arrival_gm_mean.txt", \
                "arrival_wm_wm_mean.txt"]

    pv_out_files = ["perfusion_gm_calib_masked.nii.gz", \
            "perfusion_gm_var_calib_masked.nii.gz", \
            "perfusion_wm_calib_masked.nii.gz", \
            "perfusion_wm_var_calib_masked.nii.gz", \
            "arrival_gm_masked.nii.gz", \
            "arrival_gm_var_masked.nii.gz", \
            "arrival_wm_masked.nii.gz", \
            "arrival_wm_var_masked.nii.gz", \
            "aCBV_calib.nii.gz"]

    pv_out_txt_files = ["perfusion_calib_gm_mean.txt", \
                     "perfusion_calib_wm_mean.txt", \
                     "arrival_gm_mean.txt", \
                     "arrival_wm_mean.txt"]

    surface_files = ["perfusion_calib_Atlas.dscalar.nii", \
                      "arrival_Atlas.dscalar.nii"]

    for x in nonpv_img_files:
        copy((source_path_T1 + x), (destination_path_T1 + x))
    for y in nonpv_txt_files:
        copy((source_path_T1 + y), (destination_path_T1 + y))
    for z in range(len(pv_img_files)):
        copy((source_path_T1 + pv_prefix + "/" + pv_img_files[z]), (destination_path_T1 + pv_prefix + "_" + pv_out_files[z]))
    for a in range(len(pv_txt_files)):
        copy((source_path_T1 + pv_prefix + "/" + pv_txt_files[a]), (destination_path_T1 + pv_prefix + "_" + pv_out_txt_files[a]))
    for b in surface_files:
        copy((source_path_MNI + b), (destination_path_MNI + b))
        copy((source_path_MNI_pv + b), destination_path_MNI + pv_prefix + "_" + b)
    
    cmd_cbf = ["wb_command -metric-stats", (source_path_MNI + surface_files[0]), "-reduce MEAN >", (destination_path_MNI + "perfusion_calib_cifti_mean.txt")]
   process = subprocess.Popen(cmd_cbf, stdout=subprocess.PIPE)

    cmd_pv_cbf = ["wb_command -metric-stats", (source_path_MNI_pv + surface_files[0]), "-reduce MEAN >", (destination_path_MNI + "pvcorr_perfusion_calib_cifti_mean.txt")]
    process = subprocess.Popen(cmd_pv_cbf, stdout=subprocess.PIPE)

    cmd_AAT = ["wb_command -metric-stats", (source_path_MNI + surface_files[1]), "-reduce MEAN >", (destination_path_MNI + "arrival_cifti_mean.txt")]
    process = subprocess.Popen(cmd_AAT, stdout=subprocess.PIPE)
    
    cmd_pv_AAT = ["wb_command -metric-stats", (source_path_MNI_pv + surface_files[1]), "-reduce MEAN >", (destination_path_MNI + "pvcorr_arrival_cifti_mean.txt")]
    process = subprocess.Popen(cmd_pv_AAT, stdout=subprocess.PIPE)