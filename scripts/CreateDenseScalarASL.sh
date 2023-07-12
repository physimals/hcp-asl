#!/bin/bash 
set -e
echo -e "\n START: CreateDenseScalarASL"

Subject="$1" 
AtlasDownSampleFolder="$8" #"${StudyFolder}/${Subject}/MNINonLinear/fsaverage_LR32k"  #"$1"
# T1DownSampleFolder="${StudyFolder}/${Subject}/T1w/fsaverage_LR32k"
AtlasResultsFolder="$9" #"${StudyFolder}/${Subject}/MNINonLinear/ASL/CIFTIPrepare"
T1wSpcResultsFolder="${10}" #"${StudyFolder}/${Subject}/T1w/ASL/CIFTIPrepare"
LowResMesh="$4" #"32" #"$3"
NameOfASL="$2" #perfusion_calib" #"$4"
SmoothingFWHM="$6" #"2" #"$5"
ROIFolder="$3" #"${StudyFolder}/${Subject}/MNINonLinear/ROIs" # "$6"
OutputAtlasDenseScalar="$7" #"${AtlasResultsFolder}/${NameOfASL}_Atlas" #"$7"
GrayordinatesResolution="$5" #"2" #"$8"

CARET7DIR="${11}" #"/Applications/workbench/bin_macosx64"

#Some way faster and more concise code:

${CARET7DIR}/wb_command -cifti-create-dense-scalar \
    "$OutputAtlasDenseScalar".dscalar.nii -volume \
    "$AtlasResultsFolder"/"$NameOfASL"_AtlasSubcortical_s"$SmoothingFWHM".nii.gz \
    "$ROIFolder"/Atlas_ROIs."$GrayordinatesResolution".nii.gz \
    -left-metric \
    "$T1wSpcResultsFolder"/"$NameOfASL"_s"$SmoothingFWHM".atlasroi.L."$LowResMesh"k_fs_LR.func.gii \
    -roi-left \
    "$AtlasDownSampleFolder"/"$Subject".L.atlasroi."$LowResMesh"k_fs_LR.shape.gii \
    -right-metric \
    "$T1wSpcResultsFolder"/"$NameOfASL"_s"$SmoothingFWHM".atlasroi.R."$LowResMesh"k_fs_LR.func.gii \
    -roi-right \
    "$AtlasDownSampleFolder"/"$Subject".R.atlasroi."$LowResMesh"k_fs_LR.shape.gii \

#Basic Cleanup
# rm "$NameOfASL"_AtlasSubcortical_s"$SmoothingFWHM".nii.gz
# rm "$NameOfASL"_s"$SmoothingFWHM".atlasroi.L."$LowResMesh"k_fs_LR.func.gii
# rm "$NameOfASL"_s"$SmoothingFWHM".atlasroi.R."$LowResMesh"k_fs_LR.func.gii

echo " END: CreateDenseScalarASL"