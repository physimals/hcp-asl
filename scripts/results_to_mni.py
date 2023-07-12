"""
Prepare ASL-gridded MNI-space (if needed)
Warp ASL-gridded T1w-space ASL variable to ASL-gridded MNI-space
"""

import regtricks as rt
import os.path as op
import nibabel as nb
import sys
import logging 

def main():

    logger_name = "HCPASL.results_to_mni"
    logger = logging.getLogger(logger_name)

    path_warp = sys.argv[
        1
    ]  # {StudyDir}/{SubjectID}/MNINonLinear/xfms/acpc_dc2standard.nii.gz
    path_to_T1_space_ASL_variable = sys.argv[
        2
    ]  # {StudyDir}/{SubjectID}/T1w/ASL/OxfordASL/native_space/<perfusion_variable>
    path_T1 = sys.argv[3]  # {StudyDir}/{SubjectID}/T1w/T1w_acpc_dc_restore.nii.gz
    path_MNI = sys.argv[4]  # /usr/local/fsl/data/standard/MNI152_T1_2mm.nii.gz
    path_to_lowres_MNI = sys.argv[
        5
    ]  # {StudyDir}/{SubjectID}/MNINonLinear/ASL/CIFTIPrepare/asl_grid_mni.nii.gz
    path_to_MNI_space_ASL_variable = sys.argv[
        6
    ]  # {StudyDir}/{SubjectID}/MNINonLinear/ASL/Results/<folder/perfusion_variable_MNI>

    # Make the ASL-grid MNI space target image for registration (if needed)
    if not op.isfile(path_to_lowres_MNI):
        logger.info("Creating ASL-grid MNI-space MNI image")
        perfusion_spc = rt.ImageSpace(path_to_T1_space_ASL_variable)
        mni_spc = rt.ImageSpace(path_MNI)
        mni_asl_grid = mni_spc.resize_voxels(perfusion_spc.vox_size / mni_spc.vox_size)
        nb.save(
            rt.Registration.identity().apply_to_image(path_MNI, mni_asl_grid),
            path_to_lowres_MNI,
        )
    else:
        logger.info("ASL-grid MNI-space MNI image already exists")

    # Warp ASL variable to newly prepared ASL-gridded MNI-space
    logger.info("Transforming ASL Variable to ASL-gridded MNI-space ASL")
    the_warp = rt.NonLinearRegistration.from_fnirt(path_warp, path_T1, path_MNI)
    asl_mni_mniaslgrid = the_warp.apply_to_image(
        path_to_T1_space_ASL_variable, path_to_lowres_MNI
    )
    nb.save(asl_mni_mniaslgrid, path_to_MNI_space_ASL_variable)


if __name__ == "__main__":
    main()
