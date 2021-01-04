#!/bin/bash 
set -e

################################################## OPTION PARSING #####################################################
#log_Msg "Parsing Command Line Options"

# parse arguments
Path="$1" #`opts_GetOpt1 "--path" $@`  # "$1" StudyFolder="/Users/florakennedymcconnell/Documents/Data_files/HCP/HCP_test/jack_pipeline_test" 
Subject="$2" #`opts_GetOpt1 "--subject" $@`  # "$2" SubjectID="HCA6002236"
ASLVariable="$3" #`opts_GetOpt1 "--aslvariable" $@`  # "$6" ASLVariable="perfusion_calib"
ASLVariableVar="$4"
LowResMesh="$5" #`opts_GetOpt1 "--lowresmesh" $@`  # "$6" LowResMesh=32
FinalASLResolution="$6" #`opts_GetOpt1 "--aslres" $@`  # "${14}" FinalASLResolution="2.5"
SmoothingFWHM="$7" #`opts_GetOpt1 "--smoothingFWHM" $@`  # "${14}" SmoothingFWHM="2"
GrayordinatesResolution="$8" #`opts_GetOpt1 "--grayordinatesres" $@`  # "${14}" GrayordinatesResolution="2"
RegName="$9" #`opts_GetOpt1 "--regname" $@` # RegName="MSMSulc" 
# script_path="${10}" #
CARET7DIR="${10}" #"/Users/florakennedymcconnell/Downloads/workbench/bin_macosx64" 
pvcorr="${11}"

# log_Msg "Path: ${Path}"
# log_Msg "Subject: ${Subject}"
# log_Msg "ASLVariable: ${ASLVariable}"
# log_Msg "LowResMesh: ${LowResMesh}"
# log_Msg "FinalASLResolution: ${FinalASLResolution}"
# log_Msg "SmoothingFWHM: ${SmoothingFWHM}"
# log_Msg "GrayordinatesResolution: ${GrayordinatesResolution}"
# log_Msg "RegName: ${RegName}"
# log_Msg "RUN: ${RUN}"

#Naming Conventions
AtlasSpaceFolder="MNINonLinear"
T1wFolder="T1w"
NativeFolder="Native"
ResultsFolder="Results"
DownSampleFolder="fsaverage_LR${LowResMesh}k"
ROIFolder="ROIs"
OutputAtlasDenseScalar="${ASLVariable}_Atlas"
#"/Applications/workbench/bin_macosx64"

AtlasSpaceFolder="$Path"/"$Subject"/"$AtlasSpaceFolder"
T1wFolder="$Path"/"$Subject"/"$T1wFolder"
if [ "$pvcorr" = false ] ; then
    InitialASLResults="$T1wFolder"/"ASL/TIs/OxfordASL/native_space"
    T1wSpcResultsFolder="$T1wFolder"/"ASL"/"$ResultsFolder"
    AtlasResultsFolder="$AtlasSpaceFolder"/"ASL"/"$ResultsFolder" #"$AtlasSpaceFolder"
else
    InitialASLResults="$T1wFolder"/"ASL/TIs/OxfordASL/native_space/pvcorr"
    T1wSpcResultsFolder="$T1wFolder"/"ASL"/"$ResultsFolder"/"pvcorr"
    AtlasResultsFolder="$AtlasSpaceFolder"/"ASL"/"$ResultsFolder"/"pvcorr" #"$AtlasSpaceFolder"
fi
echo "Projecting ASL Variables from: $InitialASLResults"
DownSampleFolder="$AtlasSpaceFolder"/"$DownSampleFolder"
ROIFolder="$AtlasSpaceFolder"/"$ROIFolder"

#Ribbon-based Volume to Surface mapping and resampling to standard surface

# log_Msg "Do volume to surface mapping"
# log_Msg "mkdir -p ${ResultsFolder}/OutputtoCIFTI"
mkdir -p "$AtlasResultsFolder"/OutputtoCIFTI
mkdir -p "$T1wSpcResultsFolder"/OutputtoCIFTI
VolumetoSurface.sh "$Subject" "$InitialASLResults" "$ASLVariable" \
        "$ASLVariableVar" "$T1wSpcResultsFolder"/"OutputtoCIFTI" \
        "$AtlasResultsFolder"/"OutputtoCIFTI" "$T1wFolder"/"$NativeFolder" \
        "$AtlasSpaceFolder"/"$NativeFolder" "$LowResMesh" "${RegName}" \
        "$DownSampleFolder" "$CARET7DIR"

#Surface Smoothing
# log_Msg "Surface Smoothing"
#"$script_path"/
SurfaceSmooth.sh "$Subject" "$AtlasResultsFolder"/"OutputtoCIFTI"/"$ASLVariable" \
        "$DownSampleFolder" "$LowResMesh" "$SmoothingFWHM" "$CARET7DIR"

# Transform voxelwise perfusion variables to MNI space
results_to_mni "$AtlasSpaceFolder"/"xfms"/"acpc_dc2standard.nii.gz" \
        "$InitialASLResults"/"${ASLVariable}.nii.gz" "$T1wFolder"/"T1w_acpc_dc_restore.nii.gz" \
        "/usr/local/fsl/data/standard/MNI152_T1_2mm.nii.gz" \
        "$AtlasResultsFolder"/"OutputtoCIFTI"/"asl_grid_mni.nii.gz" \
        "$AtlasResultsFolder"/"OutputtoCIFTI"/"${ASLVariable}_MNI.nii.gz"

#Subcortical Processing
# log_Msg "Subcortical Processing"
SubcorticalProcessing.sh "$ASLVariable" "$AtlasSpaceFolder" \
        "$AtlasResultsFolder"/"OutputtoCIFTI" "$FinalASLResolution" "$SmoothingFWHM" \
        "$GrayordinatesResolution" "$ROIFolder" "$CARET7DIR"

#Generation of Dense Timeseries
# log_Msg "Generation of Dense Scalar"
CreateDenseScalar.sh "$Subject" "$AtlasResultsFolder"/"OutputtoCIFTI"/"${ASLVariable}" \
        "$ROIFolder" "$LowResMesh" "$GrayordinatesResolution" "$SmoothingFWHM" \
        "$AtlasResultsFolder"/"OutputtoCIFTI"/"$OutputAtlasDenseScalar" \
        "$DownSampleFolder" "$CARET7DIR"

# log_Msg "Completed"
