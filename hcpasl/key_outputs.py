import os
from glob import glob
from pathlib import Path
from shutil import copy, move

from hcpasl.utils import sp_run


def copy_key_outputs(path, t1w_preproc, mni_raw):
    source_path_T1 = path + "/T1w/ASL/perfusion_estimation/native_space/"
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
    destination_path_MNI_voxel = (
        path + "/MNINonLinear/ASL/perfusion_estimation/std_space/"
    )

    t1_out_dir = Path(destination_path_T1)
    t1_out_dir.mkdir(exist_ok=True, parents=True)
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
        sp_run(mask_cmd)

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
        sp_run(mask_cmd)

    # Make key outputs more prominent in /T1w/ASL/ and warp the voxelwise results to
    # MNI space
    nonpv_img_files = [
        "perfusion_calib.nii.gz",
        "perfusion_var_calib.nii.gz",
        "arrival.nii.gz",
        "arrival_var.nii.gz",
        "aCBV_calib.nii.gz",
    ]
    for x in nonpv_img_files:
        copy((source_path_T1 + x), (destination_path_T1 + x))
        if x != "aCBV_calib.nii.gz":
            warp_cmd = [
                script,
                warp,
                (source_path_T1 + x),
                T1w_img,
                mni_img,
                asl_grid_mni,
                (destination_path_MNI_voxel + x),
            ]
            sp_run(warp_cmd)

    # Make key perfusion summary values more priminent in /T1w/ASL/
    nonpv_txt_files = [
        "perfusion_calib_gm_mean.txt",
        "perfusion_calib_wm_mean.txt",
        "arrival_gm_mean.txt",
        "arrival_wm_mean.txt",
    ]
    for y in nonpv_txt_files:
        copy((source_path_T1 + y), (destination_path_T1 + y))

    # Make key partial volume corrected outputs more prominent in /T1w/ASL/ and warp the
    # pvcorr voxelwise results to MNI space
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
    for en, out in zip(pv_img_files, pv_out_files):
        copy(
            (source_path_T1 + pv_prefix + "/" + en),
            (destination_path_T1 + pv_prefix + "_" + out),
        )
        if en != "aCBV_calib.nii.gz":
            pv_warp_cmd = [
                script,
                warp,
                (source_path_T1 + pv_prefix + "/" + en),
                T1w_img,
                mni_img,
                asl_grid_mni,
                (destination_path_MNI_voxel + pv_prefix + "/" + out),
            ]
            sp_run(pv_warp_cmd)

    # Make key pvcorr perfusion summary results more prominent in /T1w/ASL/
    pv_txt_files = [
        "perfusion_calib_gm_mean.txt",
        "perfusion_wm_calib_wm_mean.txt",
        "arrival_gm_mean.txt",
        "arrival_wm_wm_mean.txt",
    ]
    pv_out_txt_files = [
        "perfusion_calib_gm_mean.txt",
        "perfusion_calib_wm_mean.txt",
        "arrival_gm_mean.txt",
        "arrival_wm_mean.txt",
    ]
    for en, out in zip(pv_txt_files, pv_out_txt_files):
        copy(
            (source_path_T1 + pv_prefix + "/" + en),
            (destination_path_T1 + pv_prefix + "_" + out),
        )

    # Make key perfusion and arrival volume results more prominent in /MNINonLinear/ASL
    key_volume_non_pvec = ["perfusion_calib.nii.gz", "arrival.nii.gz"]
    key_volume_pvec = [
        "perfusion_gm_calib_masked.nii.gz",
        "perfusion_wm_calib_masked.nii.gz",
        "arrival_gm_masked.nii.gz",
        "arrival_wm_masked.nii.gz",
    ]
    for c in key_volume_non_pvec:
        copy(
            (destination_path_MNI + "perfusion_estimation/std_space/" + c),
            (destination_path_MNI + c),
        )
    for d in key_volume_pvec:
        copy(
            (
                destination_path_MNI
                + "perfusion_estimation/std_space/"
                + pv_prefix
                + "/"
                + d
            ),
            (destination_path_MNI + pv_prefix + "_" + d),
        )

    # Make key perfusion and arrival CIFTI files more prominent in /MNINonLinear/ASL/
    cifti_files = ["perfusion_calib_Atlas*.dscalar.nii", "arrival_Atlas*.dscalar.nii"]
    for b in cifti_files:
        fname = [Path(f).name for f in glob(source_path_MNI + b)]
        for f in fname:
            # Produce arrival and perfusion CIFTI summary values in /MNINonLinear/ASL/
            stem = f.replace(".dscalar.nii", "")
            cmd = [
                f"{os.environ['CARET7DIR']}/wb_command",
                "-cifti-stats",
                (source_path_MNI + f),
                "-reduce",
                "MEAN",
                ">",
                destination_path_MNI + f"{stem}_mean_nonzero.txt",
            ]
            cmd = " ".join(cmd)
            sp_run(cmd, shell=True)

            cmd = [
                f"{os.environ['CARET7DIR']}/wb_command",
                "-cifti-stats",
                (source_path_MNI_pv + f),
                "-reduce",
                "MEAN",
                ">",
                destination_path_MNI + f"pvcorr_{stem}_mean_nonzero.txt",
            ]
            cmd = " ".join(cmd)
            sp_run(cmd, shell=True)

            move((source_path_MNI + f), (destination_path_MNI + f))
            move((source_path_MNI_pv + f), destination_path_MNI + pv_prefix + "_" + f)
