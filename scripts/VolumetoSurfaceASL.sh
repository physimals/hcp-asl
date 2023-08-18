#!/bin/bash 
set -e
echo -e "\n START: RibbonVolumeToSurfaceMapping"

Subject="$1" #"${SubjectID}_V1_MR" 
ASLFolder="$2" #"$StudyFolder/$SubjectID/T1w/ASL/OxfordASL/native_space"
ASLVariable="$3" #"perfusion_calib"
ASLVariableVar="$4" # e.g. perfusion_var_calib
T1WorkingDirectory="$5" #"$StudyFolder/$SubjectID/T1w/ASL/CIFTIPrepare" #
T1wNativeFolder="$6" #"$StudyFolder/$SubjectID/T1w/Native"
AtlasSpaceNativeFolder="$7" #"$StudyFolder/$SubjectID/MNINonLinear/Native" 
LowResMesh="$8" #32
T1DownsampleFolder="$10" #"$StudyFolder/$SubjectID/T1w/fsaverage_LR32k" 
RegName="$9" #"MSMAll" 
AtlasDownSampleFolder="$11" #"$StudyFolder/$SubjectID/MNINonLinear/fsaverage_LR32k"
CARET7DIR="$12"

# Add leading underscore to regname for use in paths 
RegString=""
if [[ "$RegName" != "" ]]
then
    RegString="_$RegName"
fi
             

# Reciprocal of the ASLVariable variance is the precision
# This is required for precision-weighting of the volume to
# surface mapping.
echo " Creating ASL parameter's precision image"
${CARET7DIR}/wb_command -volume-math '1 / var' \
                            "$T1WorkingDirectory"/"$ASLVariable"_precision.nii.gz \
                            -var var \
                            "$ASLFolder"/"$ASLVariableVar".nii.gz

echo " Volume to surface mapping and resampling to atlas space"
for Hemisphere in L R ; do

  # Precision-weighted ribbon-constrained volume to surface mapping
  # of the ASL-derived variable
  ${CARET7DIR}/wb_command -volume-to-surface-mapping \
                            "$ASLFolder"/"$ASLVariable".nii.gz \
                            "$T1wNativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii \
                            "$T1WorkingDirectory"/"$ASLVariable"."$Hemisphere".native.func.gii \
                            -ribbon-constrained \
                              "$T1wNativeFolder"/"$Subject"."$Hemisphere".white.native.surf.gii \
                              "$T1wNativeFolder"/"$Subject"."$Hemisphere".pial.native.surf.gii \
                            -volume-roi \
                              "$T1WorkingDirectory"/"$ASLVariable"_precision.nii.gz \
                              -weighted \
                            -bad-vertices-out \
                              "$T1WorkingDirectory"/"$ASLVariable"."$Hemisphere".badvert_ribbonroi.native.func.gii

  # 10mm dilation
  ${CARET7DIR}/wb_command -metric-dilate \
                            "$T1WorkingDirectory"/"$ASLVariable"."$Hemisphere".native.func.gii \
                            "$T1wNativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii \
                            10 "$T1WorkingDirectory"/"$ASLVariable"."$Hemisphere".native.func.gii \
                            -bad-vertex-roi "$T1WorkingDirectory"/"$ASLVariable"."$Hemisphere".badvert_ribbonroi.native.func.gii
  
  # Mask using pre-computed ROI mask 
  ${CARET7DIR}/wb_command -metric-mask \
                           "$T1WorkingDirectory"/"$ASLVariable"."$Hemisphere".native.func.gii \
                           "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii  \
                           "$T1WorkingDirectory"/"$ASLVariable"."$Hemisphere".native.func.gii

  # Resample the surface to the 32k atlas space.
  ${CARET7DIR}/wb_command -metric-resample \
                            "$T1WorkingDirectory"/"$ASLVariable"."$Hemisphere".native.func.gii \
                            "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".sphere.${RegName}.native.surf.gii \
                            "$AtlasDownSampleFolder"/"$Subject"."$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii \
                            ADAP_BARY_AREA \
                            "$T1WorkingDirectory"/"$ASLVariable""$RegString"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.func.gii \
                            -area-surfs \
                              "$T1wNativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii \
                              "$T1DownsampleFolder"/"$Subject"."$Hemisphere".midthickness${RegString}."$LowResMesh"k_fs_LR.surf.gii \
                            -current-roi \
                              "$AtlasSpaceNativeFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii

  # 30mm dilation 
  ${CARET7DIR}/wb_command -metric-dilate \
                            "$T1WorkingDirectory"/"$ASLVariable""$RegString"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.func.gii \
                            "$T1DownsampleFolder"/"$Subject"."$Hemisphere".midthickness${RegString}."$LowResMesh"k_fs_LR.surf.gii \
                            30 "$T1WorkingDirectory"/"$ASLVariable""$RegString"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.func.gii 

  # Final masking 
  ${CARET7DIR}/wb_command -metric-mask \
                            "$T1WorkingDirectory"/"$ASLVariable""$RegString"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.func.gii \
                            "$AtlasDownSampleFolder"/"$Subject"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.shape.gii \
                            "$T1WorkingDirectory"/"$ASLVariable""$RegString"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.func.gii
done

echo " END: RibbonVolumeToSurfaceMapping"

