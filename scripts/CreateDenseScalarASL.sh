#!/bin/bash 
set -e
echo -e "\n START: CreateDenseScalarASL"

# StudyFolder="/Users/florakennedymcconnell/Documents/Data_files/HCP/HCP_test/jack_pipeline_test"
Subject="$1" #"HCA6002236" #
DownSampleFolder="$8" #"${StudyFolder}/${Subject}/MNINonLinear/fsaverage_LR32k"  #"$1"
LowResMesh="$4" #"32" #"$3"
NameOffMRI="$2" #"${StudyFolder}/${Subject}/T1w/ASL/Results/RibbonVolumetoSurface/perfusion_calib" #"$4"
SmoothingFWHM="$6" #"2" #"$5"
ROIFolder="$3" #"${StudyFolder}/${Subject}/MNINonLinear/ROIs" # "$6"
OutputAtlasDenseScalar="$7" #"${NameOffMRI}_atlas" #"$7"
GrayordinatesResolution="$5" #"2" #"$8"

CARET7DIR="$9" #"/Applications/workbench/bin_macosx64"

#Some way faster and more concise code:

${CARET7DIR}/wb_command -cifti-create-dense-scalar \
    "$OutputAtlasDenseScalar".dscalar.nii -volume \
    "$NameOffMRI"_AtlasSubcortical_s"$SmoothingFWHM".nii.gz \
    "$ROIFolder"/Atlas_ROIs."$GrayordinatesResolution".nii.gz \
    -left-metric \
    "$NameOffMRI"_s"$SmoothingFWHM".atlasroi.L."$LowResMesh"k_fs_LR.func.gii \
    -roi-left \
    "$DownSampleFolder"/"$Subject".L.atlasroi."$LowResMesh"k_fs_LR.shape.gii \
    -right-metric \
    "$NameOffMRI"_s"$SmoothingFWHM".atlasroi.R."$LowResMesh"k_fs_LR.func.gii \
    -roi-right \
    "$DownSampleFolder"/"$Subject".R.atlasroi."$LowResMesh"k_fs_LR.shape.gii \

#Basic Cleanup
rm "$NameOffMRI"_AtlasSubcortical_s"$SmoothingFWHM".nii.gz
rm "$NameOffMRI"_s"$SmoothingFWHM".atlasroi.L."$LowResMesh"k_fs_LR.func.gii
rm "$NameOffMRI"_s"$SmoothingFWHM".atlasroi.R."$LowResMesh"k_fs_LR.func.gii

echo " END: CreateDenseScalarASL"