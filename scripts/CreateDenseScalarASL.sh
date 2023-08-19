#!/bin/bash 
echo -e "\n START: CreateDenseScalarASL"

Subject="$1" 
ASLVariable="$2" #perfusion_calib"
ROIFolder="$3" #"${StudyFolder}/${Subject}/MNINonLinear/ROIs"
LowResMesh="$4" #"32"
RegName="$5" # MSMAll
GrayordinatesResolution="$6" 
SmoothingFWHM="$7" #"2"
OutputAtlasDenseScalar="$8" #"${AtlasResultsFolder}/${NameOfASL}_Atlas" 
AtlasDownSampleFolder="$9" #"${StudyFolder}/${Subject}/MNINonLinear/fsaverage_LR32k" 
AtlasResultsFolder="${10}" #"${StudyFolder}/${Subject}/MNINonLinear/ASL/CIFTIPrepare"
T1wSpcResultsFolder="${11}" #"${StudyFolder}/${Subject}/T1w/ASL/CIFTIPrepare"
CARET7DIR="${12}" #"/Applications/workbench/bin_macosx64"

# Add leading underscore to regname for use in paths 
RegString=""
if [[ "$RegName" != "" ]]
then
    RegString="_$RegName"
fi

${CARET7DIR}/wb_command -cifti-create-dense-scalar \
    "$OutputAtlasDenseScalar".dscalar.nii -volume \
    "$AtlasResultsFolder"/"$NameOfASL"_AtlasSubcortical_s"$SmoothingFWHM".nii.gz \
    "$ROIFolder"/Atlas_ROIs."$GrayordinatesResolution".nii.gz \
    -left-metric \
        "$T1wSpcResultsFolder"/"$NameOfASL""$RegString"_s"$SmoothingFWHM".atlasroi.L."$LowResMesh"k_fs_LR.func.gii \
    -roi-left \
        "$AtlasDownSampleFolder"/"$Subject".L.atlasroi."$LowResMesh"k_fs_LR.shape.gii \
    -right-metric \
        "$T1wSpcResultsFolder"/"$NameOfASL""$RegString"_s"$SmoothingFWHM".atlasroi.R."$LowResMesh"k_fs_LR.func.gii \
    -roi-right \
        "$AtlasDownSampleFolder"/"$Subject".R.atlasroi."$LowResMesh"k_fs_LR.shape.gii \

echo " END: CreateDenseScalarASL"