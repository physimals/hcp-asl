#!/bin/bash 
set -e

################################################## OPTION PARSING #####################################################
#log_Msg "Parsing Command Line Options"

# parse arguments
Path="$1" #`opts_GetOpt1 "--path" $@`  # "$1" StudyFolder="/Users/florakennedymcconnell/Documents/Data_files/HCP/HCP_test/jack_pipeline_test" 
Subject="$2" #`opts_GetOpt1 "--subject" $@`  # "$2" SubjectID="HCA6002236"
ASLVariable="$3" #`opts_GetOpt1 "--aslvariable" $@`  # "$6" ASLVariable="perfusion_calib"
LowResMesh="$4" #`opts_GetOpt1 "--lowresmesh" $@`  # "$6" LowResMesh=32
FinalASLResolution="$5" #`opts_GetOpt1 "--aslres" $@`  # "${14}" FinalASLResolution="2.5"
SmoothingFWHM="$6" #`opts_GetOpt1 "--smoothingFWHM" $@`  # "${14}" SmoothingFWHM="2"
GrayordinatesResolution="$7" #`opts_GetOpt1 "--grayordinatesres" $@`  # "${14}" GrayordinatesResolution="2"
RegName="$8" #`opts_GetOpt1 "--regname" $@` # RegName="MSMSulc" 
script_path="$9" #

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
CARET7DIR="/Applications/workbench/bin_macosx64"

AtlasSpaceFolder="$Path"/"$Subject"/"$AtlasSpaceFolder"
T1wFolder="$Path"/"$Subject"/"$T1wFolder"
T1wSpcResultsFolder="$T1wFolder"/"ASL"/"$ResultsFolder"
InitialASLResults="$T1wFolder"/"ASL/TIs/OxfordASL/native_space" 
ResultsFolder="$AtlasSpaceFolder"/"ASL"/"$ResultsFolder" #"$AtlasSpaceFolder"
DownSampleFolder="$AtlasSpaceFolder"/"$DownSampleFolder"
ROIFolder="$AtlasSpaceFolder"/"$ROIFolder"

#Ribbon-based Volume to Surface mapping and resampling to standard surface

# log_Msg "Do volume to surface mapping"
# log_Msg "mkdir -p ${ResultsFolder}/OutputtoCIFTI"
mkdir -p "$ResultsFolder"/OutputtoCIFTI
"$script_path"/VolumetoSurface.sh "$Subject" "$InitialASLResults" "$ASLVariable" \
        "$ResultsFolder"/"OutputtoCIFTI" "$T1wFolder"/"$NativeFolder" \
        "$AtlasSpaceFolder"/"$NativeFolder" "$LowResMesh" "${RegName}" \
        "$DownSampleFolder" "$CARET7DIR"

#Surface Smoothing
# log_Msg "Surface Smoothing"
"$script_path"/SurfaceSmooth.sh "$Subject" "$ResultsFolder"/"OutputtoCIFTI"/"$ASLVariable" \
        "$DownSampleFolder" "$LowResMesh" "$SmoothingFWHM" "$CARET7DIR"

# Transform voxelwise perfusion variables to MNI space
python "$script_path"/results_to_mni.py "$AtlasSpaceFolder"/"xfms"/"acpc_dc2standard.nii.gz" \
        "$InitialASLResults"/"${ASLVariable}.nii.gz" "$T1wFolder"/"T1w_acpc_dc_restore.nii.gz" \
        "/usr/local/fsl/data/standard/MNI152_T1_2mm.nii.gz" \
        "$ResultsFolder"/"OutputtoCIFTI"/"asl_grid_mni.nii.gz" \
        "$ResultsFolder"/"OutputtoCIFTI"/"${ASLVariable}_MNI.nii.gz"


#Subcortical Processing
# log_Msg "Subcortical Processing"
"$script_path"/SubcorticalProcessing.sh "$InitialASLResults" "$ASLVariable" "$AtlasSpaceFolder" \
        "$ResultsFolder"/"OutputtoCIFTI" "$FinalASLResolution" "$SmoothingFWHM" \
        "$GrayordinatesResolution" "$ROIFolder" "$CARET7DIR"

#Generation of Dense Timeseries
# log_Msg "Generation of Dense Scalar"
"$script_path"/CreateDenseScalar.sh "$Subject" "$ResultsFolder"/"OutputtoCIFTI"/"$ASLVariable" \
        "$ROIFolder" "$LowResMesh" "$GrayordinatesResolution" "$SmoothingFWHM" \
        "$ResultsFolder"/"OutputtoCIFTI"/"$OutputAtlasDenseScalar" \
        "$DownSampleFolder" "$CARET7DIR"

# log_Msg "Completed"
