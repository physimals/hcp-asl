#!/bin/bash 
set -e

################################################## OPTION PARSING #####################################################
#log_Msg "Parsing Command Line Options"

# parse arguments
Path="$1" #`opts_GetOpt1 "--path" $@`  # "$1"  
Subject="$2" #`opts_GetOpt1 "--subject" $@`  # "$2"
ASLVariable="$3" #`opts_GetOpt1 "--aslvariable" $@`  # "$6" ASLVariable="perfusion_calib"
ASLVariableVar="$4"
LowResMesh="$5" #`opts_GetOpt1 "--lowresmesh" $@`  # "$6" LowResMesh=32
FinalASLResolution="$6" #`opts_GetOpt1 "--aslres" $@`  # "${14}" FinalASLResolution="2.5"
SmoothingFWHM="$7" #`opts_GetOpt1 "--smoothingFWHM" $@`  # "${14}" SmoothingFWHM="2"
GrayordinatesResolution="$8" #`opts_GetOpt1 "--grayordinatesres" $@`  # "${14}" GrayordinatesResolution="2"
RegName="${9}" #`opts_GetOpt1 "--regname" $@` # RegName="MSMSulc" 
CARET7DIR="${10}" #"workbench/bin_macosx64 directory" 
pvcorr="${11}"
Outdir="${12}"

# log_Msg "Path: ${Path}"
# log_Msg "Subject: ${Subject}"
# log_Msg "ASLVariable: ${ASLVariable}"
# log_Msg "LowResMesh: ${LowResMesh}"
# log_Msg "FinalASLResolution: ${FinalASLResolution}"
# log_Msg "SmoothingFWHM: ${SmoothingFWHM}"
# log_Msg "GrayordinatesResolution: ${GrayordinatesResolution}"
# log_Msg "RegName: ${RegName}"
# log_Msg "RUN: ${RUN}"

#Naming 
T1wFolder="T1w"
AtlasSpaceFolder="MNINonLinear"
NativeFolder="Native"
ResultsFolder="Results"
DownSampleFolder="fsaverage_LR${LowResMesh}k"
ROIFolder="ROIs"
OutputAtlasDenseScalar="${ASLVariable}_Atlas"

AtlasSpaceFolder="$Path"/"$Subject"/"$AtlasSpaceFolder"
T1wFolder="$Path"/"$Subject"/"$T1wFolder"
ASLT1wFolder="$Path"/"$Subject"/"$Outdir"/"T1w/ASL"
if [ "$pvcorr" = false ] ; then
    InitialASLResults="$ASLT1wFolder"/"TIs/OxfordASL/native_space"
    T1wSpcResultsFolder="$Path"/"$Subject"/"$Outdir"/"T1w/ASL"/"$ResultsFolder"
    AtlasResultsFolder="$Path"/"$Subject"/"$Outdir"/"MNINonLinear/ASL"/"$ResultsFolder"
else
    InitialASLResults="$ASLT1wFolder"/"TIs/OxfordASL/native_space/pvcorr"
    T1wSpcResultsFolder="$Path"/"$Subject"/"$Outdir"/"T1w/ASL"/"$ResultsFolder"/"pvcorr"
    AtlasResultsFolder="$Path"/"$Subject"/"$Outdir"/"MNINonLinear/ASL"/"$ResultsFolder"/"pvcorr"
fi
echo "Projecting ASL Variables from: $InitialASLResults"
T1DownSampleFolder="$T1wFolder"/"$DownSampleFolder"
AtlasDownSampleFolder="$AtlasSpaceFolder"/"$DownSampleFolder"
ROIFolder="$AtlasSpaceFolder"/"$ROIFolder"

#Ribbon-based Volume to Surface mapping and resampling to standard surface

# log_Msg "Do volume to surface mapping"
# log_Msg "mkdir -p ${ResultsFolder}/OutputtoCIFTI"
mkdir -p "$AtlasResultsFolder"/OutputtoCIFTI
mkdir -p "$T1wSpcResultsFolder"/OutputtoMSMAll_32k
VolumetoSurfaceASL.sh "$Subject" "$InitialASLResults" "$ASLVariable" \
        "$ASLVariableVar" "$T1wSpcResultsFolder"/"OutputtoMSMAll_32k" \
        "$T1wFolder"/"$NativeFolder" \
        "$AtlasSpaceFolder"/"$NativeFolder" "$LowResMesh" "${RegName}" \
        "$T1DownSampleFolder" "$AtlasDownSampleFolder" "$CARET7DIR"

#Surface Smoothing
# log_Msg "Surface Smoothing"
SurfaceSmoothASL.sh "$Subject" "$T1wSpcResultsFolder"/"OutputtoMSMAll_32k"/"$ASLVariable" \
        "$T1DownSampleFolder" "$AtlasDownSampleFolder" "$LowResMesh" "$SmoothingFWHM" \
        "${RegName}" "$CARET7DIR"

# Transform voxelwise perfusion variables to MNI space
results_to_mni_asl "$AtlasSpaceFolder"/"xfms"/"acpc_dc2standard.nii.gz" \
        "$InitialASLResults"/"${ASLVariable}.nii.gz" "$T1wFolder"/"T1w_acpc_dc_restore.nii.gz" \
        "${FSLDIR}/data/standard/MNI152_T1_2mm.nii.gz" \
        "$AtlasResultsFolder"/"OutputtoCIFTI"/"asl_grid_mni.nii.gz" \
        "$AtlasResultsFolder"/"OutputtoCIFTI"/"${ASLVariable}_MNI.nii.gz"

#Subcortical Processing
# log_Msg "Subcortical Processing"
SubcorticalProcessingASL.sh "$ASLVariable" "$AtlasSpaceFolder" \
        "$AtlasResultsFolder"/"OutputtoCIFTI" "$FinalASLResolution" "$SmoothingFWHM" \
        "$GrayordinatesResolution" "$ROIFolder" "$CARET7DIR"

#Generation of Dense Timeseries
# log_Msg "Generation of Dense Scalar"
CreateDenseScalarASL.sh "$Subject" "${ASLVariable}" \
        "$ROIFolder" "$LowResMesh" "$GrayordinatesResolution" "$SmoothingFWHM" \
        "$AtlasResultsFolder"/"OutputtoCIFTI"/"$OutputAtlasDenseScalar" \
        "$AtlasDownSampleFolder" "$AtlasResultsFolder"/"OutputtoCIFTI" \
        "$T1wSpcResultsFolder"/"OutputtoMSMAll_32k" "$CARET7DIR"

# Move the T1w/ASL/Results/OutputtoCIFTI to MNINonLinear/ASL/Results/OutputtoCIFTI/T1wOutputtoCIFTI

# mkdir -p "$AtlasResultsFolder"/Native

# if [ -d "$AtlasResultsFolder"/"Native"/"T1wOutputtoCIFTI"/"OutputtoCIFTI" ]; then rm -rf "$AtlasResultsFolder"/"Native"/"T1wOutputtoCIFTI"/"OutputtoCIFTI"; fi
# mv "$T1wSpcResultsFolder"/OutputtoCIFTI "$AtlasResultsFolder"/"Native"/"T1wOutputtoCIFTI"

# log_Msg "Completed"