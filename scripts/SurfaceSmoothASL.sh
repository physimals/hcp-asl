#!/bin/bash 
set -e
echo -e "\n START: SurfaceSmoothing"

ASLVariable="$2" #$StudyFolder/$SubjectID/T1w/ASL/CIFTIPrepare/<perfusion_calib or arrival>
Subject="$1"
T1DownSampleFolder="$3" #"$StudyFolder/$SubjectID/T1w/fsaverage_LR32k"
AtlasDownSampleFolder="$4" #"$StudyFolder/$SubjectID/MNINonLinear/fsaverage_LR32k" #
LowResMesh="$5"
SmoothingFWHM="$6"
RegName="$7"

Sigma=$(echo "$SmoothingFWHM / (2 * sqrt(2 * l(2)))" | bc -l)

CARET7DIR="$8" #"workbench/bin_macosx64" directory

for Hemisphere in L R ; do

  ${CARET7DIR}/wb_command -metric-smoothing \
  "$T1DownSampleFolder"/"$Subject"."$Hemisphere".midthickness_${RegName}."$LowResMesh"k_fs_LR.surf.gii \
  "$ASLVariable"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.func.gii \
  "$Sigma" \
  "$ASLVariable"_s"$SmoothingFWHM".atlasroi."$Hemisphere"."$LowResMesh"k_fs_LR.func.gii \
  -roi \
    "$AtlasDownSampleFolder"/"$Subject"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.shape.gii

  #Basic Cleanup
  rm "$ASLVariable"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.func.gii
  
done

echo " END: SurfaceSmoothing"