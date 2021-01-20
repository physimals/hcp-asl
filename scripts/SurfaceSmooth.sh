#!/bin/bash 
set -e
echo -e "\n START: SurfaceSmoothing"

ASLVariable="$2" 
Subject="$1"
DownSampleFolder="$3" #"$StudyFolder/$SubjectID/MNINonLinear/fsaverage_LR32k" #
LowResMesh="$4"
SmoothingFWHM="$5"

Sigma=`echo "$SmoothingFWHM / ( 2 * ( sqrt ( 2 * l ( 2 ) ) ) )" | bc -l`

CARET7DIR="$6" #"workbench/bin_macosx64" directory

Subject="${Subject}_V1_MR"

for Hemisphere in L R ; do
  ${CARET7DIR}/wb_command -metric-smoothing \
  "$DownSampleFolder"/"$Subject"."$Hemisphere".midthickness."$LowResMesh"k_fs_LR.surf.gii \
  "$ASLVariable"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.func.gii \
  "$Sigma" \
  "$ASLVariable"_s"$SmoothingFWHM".atlasroi."$Hemisphere"."$LowResMesh"k_fs_LR.func.gii \
  -roi \
  "$DownSampleFolder"/"$Subject"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.shape.gii
  #Basic Cleanup
  rm "$ASLVariable"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.func.gii
done

echo " END: SurfaceSmoothing"