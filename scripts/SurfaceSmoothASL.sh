#!/bin/bash 
set -e
echo -e "\n START: SurfaceSmoothing"

Subject="$1"
ASLVariable="$2" #$StudyFolder/$SubjectID/T1w/ASL/CIFTIPrepare/<perfusion_calib or arrival>
T1DownSampleFolder="$3" #"$StudyFolder/$SubjectID/T1w/fsaverage_LR32k"
AtlasDownSampleFolder="$4" #"$StudyFolder/$SubjectID/MNINonLinear/fsaverage_LR32k" #
LowResMesh="$5"
SmoothingFWHM="$6"
RegName="$7"
CARET7DIR="$8" 

Sigma=$(echo "$SmoothingFWHM / (2 * sqrt(2 * l(2)))" | bc -l)

# Add leading underscore to regname for use in paths 
RegString=""
if [[ "$RegName" != "" ]]
then
    RegString="_$RegName"
fi

for Hemisphere in L R ; do

  ${CARET7DIR}/wb_command -metric-smoothing \
  "$T1DownSampleFolder"/"$Subject"."$Hemisphere".midthickness_${RegName}."$LowResMesh"k_fs_LR.surf.gii \
  "$ASLVariable""$RegString"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.func.gii \
  "$Sigma" \
  "$ASLVariable""$RegString"_s"$SmoothingFWHM".atlasroi."$Hemisphere"."$LowResMesh"k_fs_LR.func.gii \
  -roi \
    "$AtlasDownSampleFolder"/"$Subject"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.shape.gii

  # cleanup
  rm "$ASLVariable""$RegString"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.func.gii
  
done

echo " END: SurfaceSmoothing"