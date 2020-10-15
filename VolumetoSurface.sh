#!/bin/bash 
set -e
echo -e "\n START: RibbonVolumeToSurfaceMapping"

StudyFolder="/home/ibmeuser/HCP_test/jack_pipeline_test"
SubjectID="HCA6002236"

WorkingDirectory="$StudyFolder/$SubjectID/T1w/ASL/Results/RibbonVolumetoSurface" #"$1" # 
ASLVariable="$StudyFolder/$SubjectID/T1w/ASL/TIs/OxfordASL/native_space/perfusion_calib" #"$2" # ASL Results --> perfusion and ATT
Subject="HCA6002236_V1_MR" #"$3" # SubjectID
DownsampleFolder="$StudyFolder/$SubjectID/T1w/fsaverage_LR32k" # "$4" # 
LowResMesh=32 #"$5" #
T1wNativeFolder="$StudyFolder/$SubjectID/T1w/Native" #"$6" 
RegName=MSMSulc #"$7" # MSMSulc

CARET7DIR="/opt/workbench/bin_linux64"

mkdir -p $WorkingDirectory

NeighborhoodSmoothing="5" # May need to change
Factor="0.5" # May need to change

#Not sure what this bit is for...
LeftGreyRibbonValue="1"
RightGreyRibbonValue="1"

for Hemisphere in L R ; do
  if [ $Hemisphere = "L" ] ; then
    GreyRibbonValue="$LeftGreyRibbonValue"
  elif [ $Hemisphere = "R" ] ; then
    GreyRibbonValue="$RightGreyRibbonValue"
  fi    
  ${CARET7DIR}/wb_command -create-signed-distance-volume \
                            "$T1wNativeFolder"/"$Subject"."$Hemisphere".white.native.surf.gii \
                            "$ASLVariable".nii.gz \
                            "$WorkingDirectory"/"$Subject"."$Hemisphere".white.native.nii.gz

  ${CARET7DIR}/wb_command -create-signed-distance-volume \
                            "$T1wNativeFolder"/"$Subject"."$Hemisphere".pial.native.surf.gii \
                            "$ASLVariable".nii.gz \
                            "$WorkingDirectory"/"$Subject"."$Hemisphere".pial.native.nii.gz

  fslmaths "$WorkingDirectory"/"$Subject"."$Hemisphere".white.native.nii.gz \
                            -thr 0 -bin -mul 255 \
                            "$WorkingDirectory"/"$Subject"."$Hemisphere".white_thr0.native.nii.gz

  fslmaths "$WorkingDirectory"/"$Subject"."$Hemisphere".white_thr0.native.nii.gz \
                            -bin "$WorkingDirectory"/"$Subject"."$Hemisphere".white_thr0.native.nii.gz

  fslmaths "$WorkingDirectory"/"$Subject"."$Hemisphere".pial.native.nii.gz \
                            -uthr 0 -abs -bin -mul 255 \
                            "$WorkingDirectory"/"$Subject"."$Hemisphere".pial_uthr0.native.nii.gz

  fslmaths "$WorkingDirectory"/"$Subject"."$Hemisphere".pial_uthr0.native.nii.gz \
                            -bin "$WorkingDirectory"/"$Subject"."$Hemisphere".pial_uthr0.native.nii.gz

  fslmaths "$WorkingDirectory"/"$Subject"."$Hemisphere".pial_uthr0.native.nii.gz \
                            -mas "$WorkingDirectory"/"$Subject"."$Hemisphere".white_thr0.native.nii.gz \
                            -mul 255 "$WorkingDirectory"/"$Subject"."$Hemisphere".ribbon.nii.gz

  fslmaths "$WorkingDirectory"/"$Subject"."$Hemisphere".ribbon.nii.gz \
                            -bin -mul $GreyRibbonValue \
                            "$WorkingDirectory"/"$Subject"."$Hemisphere".ribbon.nii.gz

  rm "$WorkingDirectory"/"$Subject"."$Hemisphere".white.native.nii.gz \
    "$WorkingDirectory"/"$Subject"."$Hemisphere".white_thr0.native.nii.gz \
    "$WorkingDirectory"/"$Subject"."$Hemisphere".pial.native.nii.gz \
    "$WorkingDirectory"/"$Subject"."$Hemisphere".pial_uthr0.native.nii.gz
done

fslmaths "$WorkingDirectory"/"$Subject".L.ribbon.nii.gz \
            -add "$WorkingDirectory"/"$Subject".R.ribbon.nii.gz -bin \
            "$WorkingDirectory"/ribbon_only.nii.gz

rm "$WorkingDirectory"/"$Subject".L.ribbon.nii.gz "$WorkingDirectory"/"$Subject".R.ribbon.nii.gz

# ########## GOOD VOXELS SECTION

# fslmaths "$ASLVariable" -Tmean "$WorkingDirectory"/mean -odt float
# fslmaths "$ASLVariable" -Tstd "$WorkingDirectory"/std -odt float
# fslmaths "$WorkingDirectory"/std -div "$WorkingDirectory"/mean "$WorkingDirectory"/cov

# fslmaths "$WorkingDirectory"/cov -mas "$WorkingDirectory"/ribbon_only.nii.gz "$WorkingDirectory"/cov_ribbon

# fslmaths "$WorkingDirectory"/cov_ribbon -div `fslstats "$WorkingDirectory"/cov_ribbon -M` "$WorkingDirectory"/cov_ribbon_norm
# fslmaths "$WorkingDirectory"/cov_ribbon_norm -bin -s $NeighborhoodSmoothing "$WorkingDirectory"/SmoothNorm
# fslmaths "$WorkingDirectory"/cov_ribbon_norm -s $NeighborhoodSmoothing -div "$WorkingDirectory"/SmoothNorm -dilD "$WorkingDirectory"/cov_ribbon_norm_s$NeighborhoodSmoothing
# fslmaths "$WorkingDirectory"/cov -div `fslstats "$WorkingDirectory"/cov_ribbon -M` -div "$WorkingDirectory"/cov_ribbon_norm_s$NeighborhoodSmoothing "$WorkingDirectory"/cov_norm_modulate
# fslmaths "$WorkingDirectory"/cov_norm_modulate -mas "$WorkingDirectory"/ribbon_only.nii.gz "$WorkingDirectory"/cov_norm_modulate_ribbon

# STD=`fslstats "$WorkingDirectory"/cov_norm_modulate_ribbon -S`
# echo $STD
# MEAN=`fslstats "$WorkingDirectory"/cov_norm_modulate_ribbon -M`
# echo $MEAN
# Lower=`echo "$MEAN - ($STD * $Factor)" | bc -l`
# echo $Lower
# Upper=`echo "$MEAN + ($STD * $Factor)" | bc -l`
# echo $Upper

# fslmaths "$WorkingDirectory"/mean -bin "$WorkingDirectory"/mask
# fslmaths "$WorkingDirectory"/cov_norm_modulate -thr $Upper -bin -sub "$WorkingDirectory"/mask -mul -1 "$WorkingDirectory"/goodvoxels

# ########### END OF GOOD VOXELS SECTION

for Hemisphere in L R ; do
  # for Map in mean cov ; do
  #   ${CARET7DIR}/wb_command -volume-to-surface-mapping \
  #                           "$WorkingDirectory"/"$Map".nii.gz \
  #                           "$T1wNativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii \
  #                           "$WorkingDirectory"/"$Hemisphere"."$Map".native.func.gii \
  #                           -ribbon-constrained \
  #                           "$T1wNativeFolder"/"$Subject"."$Hemisphere".white.native.surf.gii \
  #                           "$T1wNativeFolder"/"$Subject"."$Hemisphere".pial.native.surf.gii 
  #                           # -volume-roi "$WorkingDirectory"/goodvoxels.nii.gz

  #   ${CARET7DIR}/wb_command -metric-dilate \
  #                           "$WorkingDirectory"/"$Hemisphere"."$Map".native.func.gii \
  #                           "$T1wNativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii \
  #                           10 "$WorkingDirectory"/"$Hemisphere"."$Map".native.func.gii -nearest

  #   ${CARET7DIR}/wb_command -metric-mask \
  #                           "$WorkingDirectory"/"$Hemisphere"."$Map".native.func.gii \
  #                           "$T1wNativeFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii \
  #                           "$WorkingDirectory"/"$Hemisphere"."$Map".native.func.gii

  #   ${CARET7DIR}/wb_command -volume-to-surface-mapping \
  #                           "$WorkingDirectory"/"$Map".nii.gz \
  #                           "$T1wNativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii \
  #                           "$WorkingDirectory"/"$Hemisphere"."$Map"_all.native.func.gii \
  #                           -ribbon-constrained \
  #                           "$T1wNativeFolder"/"$Subject"."$Hemisphere".white.native.surf.gii \
  #                           "$T1wNativeFolder"/"$Subject"."$Hemisphere".pial.native.surf.gii

  #   ${CARET7DIR}/wb_command -metric-mask \
  #                           "$WorkingDirectory"/"$Hemisphere"."$Map"_all.native.func.gii \
  #                           "$T1wNativeFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii \
  #                           "$WorkingDirectory"/"$Hemisphere"."$Map"_all.native.func.gii

  #   ${CARET7DIR}/wb_command -metric-resample \
  #                           "$WorkingDirectory"/"$Hemisphere"."$Map".native.func.gii \
  #                           "$T1wNativeFolder"/"$Subject"."$Hemisphere".sphere.${RegName}.native.surf.gii \
  #                           "$DownsampleFolder"/"$Subject"."$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii \
  #                           ADAP_BARY_AREA \
  #                           "$WorkingDirectory"/"$Hemisphere"."$Map"."$LowResMesh"k_fs_LR.func.gii \
  #                           -area-surfs \
  #                           "$T1wNativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii \
  #                           "$DownsampleFolder"/"$Subject"."$Hemisphere".midthickness."$LowResMesh"k_fs_LR.surf.gii \
  #                           -current-roi "$T1wNativeFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii

  #   ${CARET7DIR}/wb_command -metric-mask \
  #                           "$WorkingDirectory"/"$Hemisphere"."$Map"."$LowResMesh"k_fs_LR.func.gii \
  #                           "$DownsampleFolder"/"$Subject"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.shape.gii \
  #                           "$WorkingDirectory"/"$Hemisphere"."$Map"."$LowResMesh"k_fs_LR.func.gii

  #   ${CARET7DIR}/wb_command -metric-resample \
  #                           "$WorkingDirectory"/"$Hemisphere"."$Map"_all.native.func.gii \
  #                           "$T1wNativeFolder"/"$Subject"."$Hemisphere".sphere.${RegName}.native.surf.gii \
  #                           "$DownsampleFolder"/"$Subject"."$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii \
  #                           ADAP_BARY_AREA \
  #                           "$WorkingDirectory"/"$Hemisphere"."$Map"_all."$LowResMesh"k_fs_LR.func.gii \
  #                           -area-surfs \
  #                           "$T1wNativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii \
  #                           "$DownsampleFolder"/"$Subject"."$Hemisphere".midthickness."$LowResMesh"k_fs_LR.surf.gii \
  #                           -current-roi "$T1wNativeFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii

  #   ${CARET7DIR}/wb_command -metric-mask \
  #                           "$WorkingDirectory"/"$Hemisphere"."$Map"_all."$LowResMesh"k_fs_LR.func.gii \
  #                           "$DownsampleFolder"/"$Subject"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.shape.gii \
  #                           "$WorkingDirectory"/"$Hemisphere"."$Map"_all."$LowResMesh"k_fs_LR.func.gii
  # done

#   ${CARET7DIR}/wb_command -volume-to-surface-mapping \
#                             "$WorkingDirectory"/goodvoxels.nii.gz \
#                             "$T1wNativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii \
#                             "$WorkingDirectory"/"$Hemisphere".goodvoxels.native.func.gii \
#                             -ribbon-constrained \
#                             "$T1wNativeFolder"/"$Subject"."$Hemisphere".white.native.surf.gii \
#                             "$T1wNativeFolder"/"$Subject"."$Hemisphere".pial.native.surf.gii

#   ${CARET7DIR}/wb_command -metric-mask \
#                             "$WorkingDirectory"/"$Hemisphere".goodvoxels.native.func.gii \
#                             "$T1wNativeFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii \
#                             "$WorkingDirectory"/"$Hemisphere".goodvoxels.native.func.gii

#   ${CARET7DIR}/wb_command -metric-resample \
#                             "$WorkingDirectory"/"$Hemisphere".goodvoxels.native.func.gii \
#                             "$T1wNativeFolder"/"$Subject"."$Hemisphere".sphere.${RegName}.native.surf.gii \
#                             "$DownsampleFolder"/"$Subject"."$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii \
#                             ADAP_BARY_AREA \
#                             "$WorkingDirectory"/"$Hemisphere".goodvoxels."$LowResMesh"k_fs_LR.func.gii \
#                             -area-surfs \
#                             "$T1wNativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii \
#                             "$DownsampleFolder"/"$Subject"."$Hemisphere".midthickness."$LowResMesh"k_fs_LR.surf.gii \
#                             -current-roi \
#                             "$T1wNativeFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii

#   ${CARET7DIR}/wb_command -metric-mask \
#                             "$WorkingDirectory"/"$Hemisphere".goodvoxels."$LowResMesh"k_fs_LR.func.gii \
#                             "$DownsampleFolder"/"$Subject"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.shape.gii \
#                             "$WorkingDirectory"/"$Hemisphere".goodvoxels."$LowResMesh"k_fs_LR.func.gii

  ${CARET7DIR}/wb_command -volume-to-surface-mapping \
                            "$ASLVariable".nii.gz \
                            "$T1wNativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii \
                            "$ASLVariable"."$Hemisphere".native.func.gii \
                            -ribbon-constrained \
                            "$T1wNativeFolder"/"$Subject"."$Hemisphere".white.native.surf.gii \
                            "$T1wNativeFolder"/"$Subject"."$Hemisphere".pial.native.surf.gii 
                            # -volume-roi \
                            # "$WorkingDirectory"/goodvoxels.nii.gz

  ${CARET7DIR}/wb_command -metric-dilate \
                            "$ASLVariable"."$Hemisphere".native.func.gii \
                            "$T1wNativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii \
                            10 "$ASLVariable"."$Hemisphere".native.func.gii -nearest

  ${CARET7DIR}/wb_command -metric-mask \
                            "$ASLVariable"."$Hemisphere".native.func.gii \
                            "$T1wNativeFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii  \
                            "$ASLVariable"."$Hemisphere".native.func.gii

  ${CARET7DIR}/wb_command -metric-resample \
                            "$ASLVariable"."$Hemisphere".native.func.gii \
                            "$T1wNativeFolder"/"$Subject"."$Hemisphere".sphere.${RegName}.native.surf.gii \
                            "$DownsampleFolder"/"$Subject"."$Hemisphere".sphere."$LowResMesh"k_fs_LR.surf.gii \
                            ADAP_BARY_AREA "$ASLVariable"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.func.gii \
                            -area-surfs \
                            "$T1wNativeFolder"/"$Subject"."$Hemisphere".midthickness.native.surf.gii \
                            "$DownsampleFolder"/"$Subject"."$Hemisphere".midthickness."$LowResMesh"k_fs_LR.surf.gii \
                            -current-roi \
                            "$T1wNativeFolder"/"$Subject"."$Hemisphere".roi.native.shape.gii

  ${CARET7DIR}/wb_command -metric-dilate \
                            "$ASLVariable"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.func.gii \
                            "$DownsampleFolder"/"$Subject"."$Hemisphere".midthickness."$LowResMesh"k_fs_LR.surf.gii \
                            30 "$ASLVariable"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.func.gii -nearest

  ${CARET7DIR}/wb_command -metric-mask \
                            "$ASLVariable"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.func.gii \
                            "$DownsampleFolder"/"$Subject"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.shape.gii \
                            "$ASLVariable"."$Hemisphere".atlasroi."$LowResMesh"k_fs_LR.func.gii
done

echo " END: RibbonVolumeToSurfaceMapping"

