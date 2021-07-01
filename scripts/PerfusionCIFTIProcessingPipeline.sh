#!/bin/bash 
set -e

################################################## OPTION PARSING #####################################################
#log_Msg "Parsing Command Line Options"

# parse arguments
Path="$1" #`opts_GetOpt1 "--path" $@`  # "$1"  
Subject="$2" #`opts_GetOpt1 "--subject" $@`  # "$2" 
Visit="$3"
ASLVariable="$4" #`opts_GetOpt1 "--aslvariable" $@`  # "$6" ASLVariable="perfusion_calib"
ASLVariableVar="$5"
LowResMesh="$6" #`opts_GetOpt1 "--lowresmesh" $@`  # "$6" LowResMesh=32
FinalASLResolution="$7" #`opts_GetOpt1 "--aslres" $@`  # "${14}" FinalASLResolution="2.5"
SmoothingFWHM="$8" #`opts_GetOpt1 "--smoothingFWHM" $@`  # "${14}" SmoothingFWHM="2"
GrayordinatesResolution="$9" #`opts_GetOpt1 "--grayordinatesres" $@`  # "${14}" GrayordinatesResolution="2"
RegName="${10}" #`opts_GetOpt1 "--regname" $@` # RegName="MSMSulc" 
CARET7DIR="${11}" #"workbench/bin_macosx64 directory" 
pvcorr="${12}"
Outdir="${13}"

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
VisitDir="${Subject}_${Visit}_MR"
T1wFolder="${VisitDir}/T1w"
AtlasSpaceFolder="${VisitDir}/MNINonLinear"
NativeFolder="Native"
ResultsFolder="Results"
DownSampleFolder="fsaverage_LR${LowResMesh}k"
ROIFolder="ROIs"
OutputAtlasDenseScalar="${ASLVariable}_Atlas"

AtlasSpaceFolder="$Path"/"$Subject"/"$AtlasSpaceFolder"
T1wFolder="$Path"/"$Subject"/"$T1wFolder"
ASLT1wFolder="$Path"/"$Subject"/"$VisitDir"/"$Outdir"/"T1w/ASL"
if [ "$pvcorr" = false ] ; then
    InitialASLResults="$ASLT1wFolder"/"TIs/OxfordASL/native_space"
    T1wSpcResultsFolder="$Path"/"$Subject"/"$VisitDir"/"$Outdir"/"T1w/ASL"/"$ResultsFolder"
    AtlasResultsFolder="$Path"/"$Subject"/"$VisitDir"/"$Outdir"/"MNINonLinear/ASL"/"$ResultsFolder"
else
    InitialASLResults="$ASLT1wFolder"/"TIs/OxfordASL/native_space/pvcorr"
    T1wSpcResultsFolder="$Path"/"$Subject"/"$VisitDir"/"$Outdir"/"T1w/ASL"/"$ResultsFolder"/"pvcorr"
    AtlasResultsFolder="$Path"/"$Subject"/"$VisitDir"/"$Outdir"/"MNINonLinear/ASL"/"$ResultsFolder"/"pvcorr"
fi
echo "Projecting ASL Variables from: $InitialASLResults"
DownSampleFolder="$AtlasSpaceFolder"/"$DownSampleFolder"
ROIFolder="$AtlasSpaceFolder"/"$ROIFolder"

#Ribbon-based Volume to Surface mapping and resampling to standard surface

# log_Msg "Do volume to surface mapping"
# log_Msg "mkdir -p ${ResultsFolder}/OutputtoCIFTI"
mkdir -p "$AtlasResultsFolder"/OutputtoCIFTI
mkdir -p "$T1wSpcResultsFolder"/OutputtoCIFTI
VolumetoSurface.sh "$VisitDir" "$InitialASLResults" "$ASLVariable" \
        "$ASLVariableVar" "$T1wSpcResultsFolder"/"OutputtoCIFTI" \
        "$AtlasResultsFolder"/"OutputtoCIFTI" "$T1wFolder"/"$NativeFolder" \
        "$AtlasSpaceFolder"/"$NativeFolder" "$LowResMesh" "${RegName}" \
        "$DownSampleFolder" "$CARET7DIR"

#Surface Smoothing
# log_Msg "Surface Smoothing"
SurfaceSmooth.sh "$VisitDir" "$AtlasResultsFolder"/"OutputtoCIFTI"/"$ASLVariable" \
        "$DownSampleFolder" "$LowResMesh" "$SmoothingFWHM" "$CARET7DIR"

# Transform voxelwise perfusion variables to MNI space
results_to_mni "$AtlasSpaceFolder"/"xfms"/"acpc_dc2standard.nii.gz" \
        "$InitialASLResults"/"${ASLVariable}.nii.gz" "$T1wFolder"/"T1w_acpc_dc_restore.nii.gz" \
        "${FSLDIR}/data/standard/MNI152_T1_2mm.nii.gz" \
        "$AtlasResultsFolder"/"OutputtoCIFTI"/"asl_grid_mni.nii.gz" \
        "$AtlasResultsFolder"/"OutputtoCIFTI"/"${ASLVariable}_MNI.nii.gz"

#Subcortical Processing
# log_Msg "Subcortical Processing"
SubcorticalProcessing.sh "$ASLVariable" "$AtlasSpaceFolder" \
        "$AtlasResultsFolder"/"OutputtoCIFTI" "$FinalASLResolution" "$SmoothingFWHM" \
        "$GrayordinatesResolution" "$ROIFolder" "$CARET7DIR"

#Generation of Dense Timeseries
# log_Msg "Generation of Dense Scalar"
CreateDenseScalar.sh "$VisitDir" "$AtlasResultsFolder"/"OutputtoCIFTI"/"${ASLVariable}" \
        "$ROIFolder" "$LowResMesh" "$GrayordinatesResolution" "$SmoothingFWHM" \
        "$AtlasResultsFolder"/"OutputtoCIFTI"/"$OutputAtlasDenseScalar" \
        "$DownSampleFolder" "$CARET7DIR"

# Move the T1w/ASL/Results/OutputtoCIFTI to MNINonLinear/ASL/Results/OutputtoCIFTI/T1wOutputtoCIFTI

mkdir -p "$AtlasResultsFolder"/Native
mv "$T1wSpcResultsFolder"/OutputtoCIFTI "$AtlasResultsFolder"/"Native"/"T1wOutputtoCIFTI"

# log_Msg "Completed"
