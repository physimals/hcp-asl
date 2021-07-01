#!/bin/bash 
set -e
script_name="SubcorticalProcessing.sh"
echo "${script_name}: START"

AtlasSpaceFolder="$2" #"${StudyFolder}/${SubjectID}/${SubjectID}_${Visit}_MR/MNINonLinear" # 
echo "${script_name}: AtlasSpaceFolder: ${AtlasSpaceFolder}"

ROIFolder="$7" #"${AtlasSpaceFolder}/ROIs" # 
echo "${script_name}: ROIFolder: ${ROIFolder}"

FinalASLResolution="$4" # "2.5" # 
echo "${script_name}: FinalASLResolution: ${FinalASLResolution}"

ResultsFolder="$3" #"${StudyFolder}/${SubjectID}/${SubjectID}_${Visit}_MR/T1w/ASL/Results/RibbonVolumetoSurface" # 
echo "${script_name}: ResultsFolder: ${ResultsFolder}"

ASLVariable="$1" #"perfusion_calib" #
echo "${script_name}: ASLVariable: ${ASLVariable}"

Visit="$9" # visit, i.e. V1 or V2
echo "${script_name}: Visit: ${Visit}"

SmoothingFWHM="$5" #"2" #
echo "${script_name}: SmoothingFWHM: ${SmoothingFWHM}"

BrainOrdinatesResolution="$6" #"2" #
echo "${script_name}: BrainOrdinatesResolution: ${BrainOrdinatesResolution}"

VolumeASL="${ResultsFolder}/${ASLVariable}"
echo "${script_name}: VolumeASL: ${VolumeASL}"

Sigma=`echo "$SmoothingFWHM / ( 2 * ( sqrt ( 2 * l ( 2 ) ) ) )" | bc -l`
echo "${script_name}: Sigma: ${Sigma}"

CARET7DIR="$8" #"/Applications/workbench/bin_macosx64"
#NOTE: wmparc has dashes in structure names, which -cifti-create-* won't accept
#ROIs files have acceptable structure names

#deal with fsl_sub being silly when we want to use numeric equality on decimals
unset POSIXLY_CORRECT

#### NEED TO TRANSFORM SUBCORTICAL METRICS TO MNI IN VOLUME SPACE
#### This can be done using results_to_mni.py

#generate subject-roi space fMRI cifti for subcortical
if [[ `echo "$BrainOrdinatesResolution == $FinalASLResolution" | bc -l | cut -f1 -d.` == "1" ]]
then
    echo "${script_name}: Creating subject-roi subcortical cifti at same resolution as output"
    ${CARET7DIR}/wb_command -cifti-create-dense-scalar \
        ${ResultsFolder}/${ASLVariable}_temp_subject.dscalar.nii \
        -volume "${VolumeASL}_MNI".nii.gz \
        "$ROIFolder"/ROIs."$BrainOrdinatesResolution".nii.gz
else
    echo "${script_name}: Creating subject-roi subcortical cifti at differing fMRI resolution"
    ${CARET7DIR}/wb_command -volume-affine-resample \
        "$ROIFolder"/ROIs."$BrainOrdinatesResolution".nii.gz \
        $FSLDIR/etc/flirtsch/ident.mat \
        "${VolumeASL}_MNI".nii.gz ENCLOSING_VOXEL \
        "$ResultsFolder"/ROIs."$FinalASLResolution".nii.gz

    ${CARET7DIR}/wb_command -cifti-create-dense-scalar \
        ${ResultsFolder}/${ASLVariable}_temp_subject.dscalar.nii \
        -volume "${VolumeASL}_MNI".nii.gz \
        "$ResultsFolder"/ROIs."$FinalASLResolution".nii.gz

    rm -f "$ResultsFolder"/ROIs."$FinalASLResolution".nii.gz
fi

echo "${script_name}: Dilating out zeros"
#dilate out any exact zeros in the input data, for instance if the brain mask is wrong
${CARET7DIR}/wb_command -cifti-dilate \
    ${ResultsFolder}/${ASLVariable}_temp_subject.dscalar.nii \
    COLUMN 0 10 \
    ${ResultsFolder}/${ASLVariable}_temp_subject_dilate.dscalar.nii

rm -f ${ResultsFolder}/${ASLVariable}_temp_subject.dscalar.nii

echo "${script_name}: Generate atlas subcortical template cifti"
${CARET7DIR}/wb_command -cifti-create-label \
    ${ResultsFolder}/${ASLVariable}_temp_template.dlabel.nii \
    -volume "$ROIFolder"/Atlas_ROIs."$BrainOrdinatesResolution".nii.gz \
    "$ROIFolder"/Atlas_ROIs."$BrainOrdinatesResolution".nii.gz

if [[ `echo "${Sigma} > 0" | bc -l | cut -f1 -d.` == "1" ]]
then
    echo "${script_name}: Smoothing and resampling"
    #this is the whole timeseries, so don't overwrite, in order to allow on-disk writing, then delete temporary
    ${CARET7DIR}/wb_command -cifti-smoothing \
        ${ResultsFolder}/${ASLVariable}_temp_subject_dilate.dscalar.nii 0 ${Sigma} \
        COLUMN ${ResultsFolder}/${ASLVariable}_temp_subject_smooth.dscalar.nii \
        -fix-zeros-volume
    #resample, delete temporary
    ${CARET7DIR}/wb_command -cifti-resample \
        ${ResultsFolder}/${ASLVariable}_temp_subject_smooth.dscalar.nii \
        COLUMN ${ResultsFolder}/${ASLVariable}_temp_template.dlabel.nii \
        COLUMN ADAP_BARY_AREA CUBIC ${ResultsFolder}/${ASLVariable}_temp_atlas.dscalar.nii \
        -volume-predilate 10

    rm -f ${ResultsFolder}/${ASLVariable}_temp_subject_smooth.dscalar.nii
else
    echo "${script_name}: Resampling"
    ${CARET7DIR}/wb_command -cifti-resample \
        ${ResultsFolder}/${ASLVariable}_temp_subject_dilate.dscalar.nii \
        COLUMN ${ResultsFolder}/${ASLVariable}_temp_template.dlabel.nii \
        COLUMN ADAP_BARY_AREA CUBIC ${ResultsFolder}/${ASLVariable}_temp_atlas.dscalar.nii \
        -volume-predilate 10
fi

#delete common temporaries
rm -f ${ResultsFolder}/${ASLVariable}_temp_subject_dilate.dscalar.nii
rm -f ${ResultsFolder}/${ASLVariable}_temp_template.dlabel.nii

#write output volume, delete temporary
#NOTE: $VolumeASL contains a path in it, it is not a file in the current directory
${CARET7DIR}/wb_command -cifti-separate \
    ${ResultsFolder}/${ASLVariable}_temp_atlas.dscalar.nii COLUMN \
    -volume-all ${ResultsFolder}/${ASLVariable}_AtlasSubcortical_s"$SmoothingFWHM".nii.gz

rm -f ${ResultsFolder}/${ASLVariable}_temp_atlas.dscalar.nii

echo "${script_name}: END"

