import os
import subprocess
from pathlib import Path
from shutil import copy


def copy_key_outputs(path, t1w_preproc, mni_raw):
    source_path_T1 = path + "/T1w/ASL/OxfordASL/native_space/"
    destination_path_T1 = path + "/T1w/ASL/"

    source_path_MNI = path + "/MNINonLinear/ASL/CIFTIPrepare/"
    source_path_MNI_pv = path + "/MNINonLinear/ASL/CIFTIPrepare/pvcorr/"
    destination_path_MNI = path + "/MNINonLinear/ASL/"

    pv_prefix = "pvcorr"

    # Set up path variables needed to warp voxelwise results to MNI space
    script = "results_to_mni_asl"
    warp = mni_raw + "/xfms/acpc_dc2standard.nii.gz"
    T1w_img = t1w_preproc + "/T1w_acpc_dc_restore.nii.gz"
    fsldir = Path(os.environ["FSLDIR"])
    mni_img = fsldir / "data/standard/MNI152_T1_2mm.nii.gz"
    asl_grid_mni = path + "/MNINonLinear/ASL/CIFTIPrepare/asl_grid_mni.nii.gz"
    destination_path_MNI_voxel = path + "/MNINonLinear/ASL/OxfordASL/std_space/"

    t1_out_dir = Path(destination_path_T1)
    t1_out_dir.mkdir(exist_ok=True)
    out_MNI_dir = Path(destination_path_MNI)
    out_MNI_dir.mkdir(exist_ok=True, parents=True)
    out_MNI_dir_voxel_pvcorr = Path(destination_path_MNI_voxel) / pv_prefix
    out_MNI_dir_voxel_pvcorr.mkdir(exist_ok=True, parents=True)

    # mask pvcorr parameter variance estimates
    gm_mask = source_path_T1 + "gm_mask.nii.gz"
    gm_pvcorr_vars = ["perfusion_var_calib.nii.gz", "arrival_var.nii.gz"]
    gm_pvcorr_vars_out = [
        "perfusion_var_calib_masked.nii.gz",
        "arrival_var_masked.nii.gz",
    ]

    # Mask grey matter partial volume corrected perfusion and arrival results to restrict
    # to grey matter-contaiming voxels only
    for gm_pvcorr_var, gm_pvcorr_var_out in zip(gm_pvcorr_vars, gm_pvcorr_vars_out):
        mask_cmd = [
            "fslmaths",
            source_path_T1 + pv_prefix + "/" + gm_pvcorr_var,
            "-mas",
            gm_mask,
            source_path_T1 + pv_prefix + "/" + gm_pvcorr_var_out,
        ]
        process = subprocess.Popen(mask_cmd, stdout=subprocess.PIPE)

    wm_mask = source_path_T1 + "wm_mask.nii.gz"
    wm_pvcorr_vars = ["perfusion_wm_var_calib.nii.gz", "arrival_wm_var.nii.gz"]
    wm_pvcorr_vars_out = [
        "perfusion_wm_var_calib_masked.nii.gz",
        "arrival_wm_var_masked.nii.gz",
    ]

    # Mask white matter partial volume corrected perfusion and arrival results to restrict
    # to white matter-containing voxels only
    for wm_pvcorr_var, wm_pvcorr_var_out in zip(wm_pvcorr_vars, wm_pvcorr_vars_out):
        mask_cmd = [
            "fslmaths",
            source_path_T1 + pv_prefix + "/" + wm_pvcorr_var,
            "-mas",
            wm_mask,
            source_path_T1 + pv_prefix + "/" + wm_pvcorr_var_out,
        ]
        process = subprocess.Popen(mask_cmd, stdout=subprocess.PIPE)

    nonpv_img_files = [
        "perfusion_calib.nii.gz",
        "perfusion_var_calib.nii.gz",
        "arrival.nii.gz",
        "arrival_var.nii.gz",
        "aCBV_calib.nii.gz",
    ]

    nonpv_txt_files = [
        "perfusion_calib_gm_mean.txt",
        "perfusion_calib_wm_mean.txt",
        "arrival_gm_mean.txt",
        "arrival_wm_mean.txt",
    ]

    pv_img_files = [
        "perfusion_calib_masked.nii.gz",
        "perfusion_var_calib_masked.nii.gz",
        "perfusion_wm_calib_masked.nii.gz",
        "perfusion_wm_var_calib_masked.nii.gz",
        "arrival_masked.nii.gz",
        "arrival_var_masked.nii.gz",
        "arrival_wm_masked.nii.gz",
        "arrival_wm_var_masked.nii.gz",
        "aCBV_calib.nii.gz",
    ]

    pv_txt_files = [
        "perfusion_calib_gm_mean.txt",
        "perfusion_wm_calib_wm_mean.txt",
        "arrival_gm_mean.txt",
        "arrival_wm_wm_mean.txt",
    ]

    # Rename GM pvcorr results to clarify that they are grey matter-related
    pv_out_files = [
        "perfusion_gm_calib_masked.nii.gz",
        "perfusion_gm_var_calib_masked.nii.gz",
        "perfusion_wm_calib_masked.nii.gz",
        "perfusion_wm_var_calib_masked.nii.gz",
        "arrival_gm_masked.nii.gz",
        "arrival_gm_var_masked.nii.gz",
        "arrival_wm_masked.nii.gz",
        "arrival_wm_var_masked.nii.gz",
        "aCBV_calib.nii.gz",
    ]

    pv_out_txt_files = [
        "perfusion_calib_gm_mean.txt",
        "perfusion_calib_wm_mean.txt",
        "arrival_gm_mean.txt",
        "arrival_wm_mean.txt",
    ]

    surface_files = ["perfusion_calib_Atlas.dscalar.nii", "arrival_Atlas.dscalar.nii"]

    key_volume_non_pvec = ["perfusion_calib.nii.gz", "arrival.nii.gz"]

    key_volume_pvec = [
        "perfusion_gm_calib_masked.nii.gz",
        "perfusion_wm_calib_masked.nii.gz",
        "arrival_gm_masked.nii.gz",
        "arrival_wm_masked.nii.gz",
    ]

    # Make key outputs more prominent in /T1w/ASL/ and warp the voxelwise results to
    # MNI space
    for x in nonpv_img_files:
        copy((source_path_T1 + x), (destination_path_T1 + x))
        if x is not "aCBV_calib.nii.gz":
            warp_cmd = [
                script,
                warp,
                (source_path_T1 + x),
                T1w_img,
                mni_img,
                asl_grid_mni,
                (destination_path_MNI_voxel + x),
            ]
            subprocess.run(warp_cmd, check=True)

    # Make key perfusion summary values more priminent in /T1w/ASL/
    for y in nonpv_txt_files:
        copy((source_path_T1 + y), (destination_path_T1 + y))

    # Make key partial volume corrected outputs more prominent in /T1w/ASL/ and warp the
    # pvcorr voxelwise results to MNI space
    for z in range(len(pv_img_files)):
        copy(
            (source_path_T1 + pv_prefix + "/" + pv_img_files[z]),
            (destination_path_T1 + pv_prefix + "_" + pv_out_files[z]),
        )
        if pv_img_files[z] is not "aCBV_calib.nii.gz":
            pv_warp_cmd = [
                script,
                warp,
                (source_path_T1 + pv_prefix + "/" + pv_img_files[z]),
                T1w_img,
                mni_img,
                asl_grid_mni,
                (destination_path_MNI_voxel + pv_prefix + "/" + pv_out_files[z]),
            ]
            subprocess.run(pv_warp_cmd, check=True)

    # Make key pvcorr perfusion summary results more prominent in /T1w/ASL/
    for a in range(len(pv_txt_files)):
        copy(
            (source_path_T1 + pv_prefix + "/" + pv_txt_files[a]),
            (destination_path_T1 + pv_prefix + "_" + pv_out_txt_files[a]),
        )

    # Make key perfusion and arrival CIFTI files more prominent in /MNINonLinear/ASL/
    for b in surface_files:
        copy((source_path_MNI + b), (destination_path_MNI + b))
        copy((source_path_MNI_pv + b), destination_path_MNI + pv_prefix + "_" + b)

    # Make key perfusion and arrival volume results more prominent in /MNINonLinear/ASL
    for c in key_volume_non_pvec:
        copy(
            (destination_path_MNI + "OxfordASL/std_space/" + c),
            (destination_path_MNI + c),
        )
    for d in key_volume_pvec:
        copy(
            (destination_path_MNI + "OxfordASL/std_space/" + pv_prefix + "/" + d),
            (destination_path_MNI + pv_prefix + "_" + d),
        )

    # Produce arrival and perfusion CIFTI summary values in /MNINonLinear/ASL/
    cmd_cbf = [
        "wb_command",
        "-cifti-stats",
        (source_path_MNI + surface_files[0]),
        "-reduce",
        "MEAN",
        ">",
        destination_path_MNI + "perfusion_calib_cifti_mean_nonzero.txt",
    ]
    subprocess.run(" ".join(cmd_cbf), shell=True)

    cmd_pv_cbf = [
        "wb_command",
        "-cifti-stats",
        (source_path_MNI_pv + surface_files[0]),
        "-reduce",
        "MEAN",
        ">",
        destination_path_MNI + "pvcorr_perfusion_calib_cifti_mean_nonzero.txt",
    ]
    subprocess.run(" ".join(cmd_pv_cbf), shell=True)

    cmd_AAT = [
        "wb_command",
        "-cifti-stats",
        (source_path_MNI + surface_files[1]),
        "-reduce",
        "MEAN",
        ">",
        destination_path_MNI + "arrival_cifti_mean_nonzero.txt",
    ]
    subprocess.run(" ".join(cmd_AAT), shell=True)

    cmd_pv_AAT = [
        "wb_command",
        "-cifti-stats",
        (source_path_MNI_pv + surface_files[1]),
        "-reduce",
        "MEAN",
        ">",
        destination_path_MNI + "pvcorr_arrival_cifti_mean_nonzero.txt",
    ]
    subprocess.run(" ".join(cmd_pv_AAT), shell=True)
