#!/bin/bash 
set -e

# parse arguments
Path="$1"
Subject="$2"
ASLVariable="$3" 
ASLVariableVar="$4"
LowResMesh="$5" # 32
FinalASLResolution="$6" #2.5
SmoothingFWHM="$7" # "2"
GrayordinatesResolution="$8" # "2"
RegName="$9" # "MSMAll" 
CARET7DIR="${10}" #"workbench/bin_macosx64 directory" 
pvcorr="${11}"
Outdir="${12}"

#Naming 
T1wFolder="T1w"
AtlasSpaceFolder="MNINonLinear"
NativeFolder="Native"
ResultsFolder="CIFTIPrepare"
DownSampleFolder="fsaverage_LR${LowResMesh}k"
ROIFolder="ROIs"
RegString=""
if [[ "$RegName" != "" ]]
then
    RegString="_$RegName"
fi

AtlasSpaceFolder="$Path"/"$Subject"/"$AtlasSpaceFolder"
T1wFolder="$Path"/"$Subject"/"$T1wFolder"
ASLT1wFolder="$Path"/"$Subject"/"$Outdir"/"T1w/ASL"
if [ "$pvcorr" = false ] ; then
    InitialASLResults="$ASLT1wFolder"/"OxfordASL/native_space"
    T1wSpcResultsFolder="$Path"/"$Subject"/"$Outdir"/"T1w/ASL"/"$ResultsFolder"
    AtlasResultsFolder="$Path"/"$Subject"/"$Outdir"/"MNINonLinear/ASL"/"$ResultsFolder"
else
    InitialASLResults="$ASLT1wFolder"/"OxfordASL/native_space/pvcorr"
    T1wSpcResultsFolder="$Path"/"$Subject"/"$Outdir"/"T1w/ASL"/"$ResultsFolder"/"pvcorr"
    AtlasResultsFolder="$Path"/"$Subject"/"$Outdir"/"MNINonLinear/ASL"/"$ResultsFolder"/"pvcorr"
fi
echo "Projecting ASL Variables from: $InitialASLResults"
T1DownSampleFolder="$T1wFolder"/"$DownSampleFolder"
AtlasDownSampleFolder="$AtlasSpaceFolder"/"$DownSampleFolder"
ROIFolder="$AtlasSpaceFolder"/"$ROIFolder"

#Ribbon-based Volume to Surface mapping and resampling to standard surface

# log_Msg "Do volume to surface mapping"
mkdir -p "$T1wSpcResultsFolder"
VolumetoSurfaceASL.sh "$Subject" "$InitialASLResults" "$ASLVariable" \
        "$ASLVariableVar" "$T1wSpcResultsFolder" \
        "$T1wFolder"/"$NativeFolder" \
        "$AtlasSpaceFolder"/"$NativeFolder" "$LowResMesh" "$RegName" \
        "$T1DownSampleFolder" "$AtlasDownSampleFolder" "$CARET7DIR"

#Surface Smoothing
SurfaceSmoothASL.sh "$Subject" "$T1wSpcResultsFolder"/"$ASLVariable" \
        "$T1DownSampleFolder" "$AtlasDownSampleFolder" "$LowResMesh" "$SmoothingFWHM" \
        "$RegName" "$CARET7DIR"

# Transform voxelwise perfusion variables to MNI space
mkdir -p "$AtlasResultsFolder"
results_to_mni_asl "$AtlasSpaceFolder"/"xfms"/"acpc_dc2standard.nii.gz" \
        "$InitialASLResults"/"${ASLVariable}.nii.gz" "$T1wFolder"/"T1w_acpc_dc_restore.nii.gz" \
        "${FSLDIR}/data/standard/MNI152_T1_2mm.nii.gz" \
        "$AtlasResultsFolder"/"asl_grid_mni.nii.gz" \
        "$AtlasResultsFolder"/"${ASLVariable}_MNI.nii.gz"

#Subcortical Processing
SubcorticalProcessingASL.sh "$ASLVariable" "$AtlasSpaceFolder" \
        "$AtlasResultsFolder" "$FinalASLResolution" "$SmoothingFWHM" \
        "$GrayordinatesResolution" "$ROIFolder" "$CARET7DIR"

#Generation of dense scalar image 
OutputAtlasDenseScalar="${ASLVariable}_Atlas${RegString}"
CreateDenseScalarASL.sh "$Subject" "$ASLVariable" \
        "$ROIFolder" "$LowResMesh" "$RegName" "$GrayordinatesResolution" \
        "$SmoothingFWHM" "$AtlasResultsFolder"/"$OutputAtlasDenseScalar" \
        "$AtlasDownSampleFolder" "$AtlasResultsFolder" \
        "$T1wSpcResultsFolder" "$CARET7DIR"
