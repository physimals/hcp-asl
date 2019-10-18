
# Simple ASL pipeline for the Human Connectome Project LifeSpan Aging 
# extension

# This pipeline is based on standard volume-based ASL perfusion analysis
# followed by projection of the resulting perfusion estimates 
# onto the cortical surface where appropriate i.e. generation of a CIFTI
# representation of the perfusion results.

# One feature of the HCP minimal processing pipeline is gradient 
# non-linearity correction. This is not standardly carried out in the 
# oxford_asl analysis tools, and so will be incorporated here prior to 
# BASIL analysis.


######################### Analysis Steps ###############################

### STEP 1
# Standard structural analysis pipeline: currently FSL_ANAT but will 
# gradually move towards using outputs from the HCP structural pipelines

fsl_anat -i <input_data> -o <anat_dir>

### STEP 2
# Run gradient distortion correction (gradunwarp) to generate GDC warp

fslroi <asl_data> <asl_data_vol_1> 0 1

InputCoeffs="/home/ibmeuser/projects/Pipelines/global/config/coeff_AS82_Prisma.grad"

# be careful about paths and working directory here - gradient_unwarp will
# output files directly to the working directory
####### !!!!! Also there is a python version conflict here !!!!! #######
gradient_unwarp.py <asl_data_vol_1> <grad_out> siemens -g $InputCoeffs -n

convertwarp --abs --ref=<grad_out> --warp1=fullWarp_abs.nii.gz --relout --out=<gdc_warp>

### STEP 3
# Run topup to obtain fieldmap required to correct for susceptibility
# induced distortions

fsl_merge -t <merged_fieldmaps> <PA_fieldmap> <AP_fieldmap>

mcflirt -in <merged_fieldmaps> 

echo "0 1 0 0.04845" > <topup_params.txt>
echo "0 -1 0 0.04845" >> <topup_params.txt>

topup --imain=<merged_fieldmaps_mcf> --datain=<topup_params.txt> \
    --config=b02b0.cnf --out=<topupresult> --fout=<topupresult_fmap>

# Average then brain extract distortion corrected fieldmaps for input to oxford_asl

fslmaths <topupresult> -Tmean fmapmag

bet fmapmag fmapmagbrain

# CHECK FIELDMAP PROCESSING 
# Fieldmap provided from topup is in Hz, and must be converted to units
# of rads/s

fslmaths <topupresult_fmap> -mul 3.14159 -mul 2 <fmap>

### STEP 4
# Run oxford_asl for perfusion analysis, including the outputs from
# each of the previous steps

fslroi <asl_data_file> <asl_data> 0 <end_of_data>

fslroi <asl_data_file> <calibration_images> <end_point>

oxford_asl -i <asl_data> \
    -o <output_filename> \
    --casl \
    --iaf=tc \
    --ibf=tis \
    --tis=1.7,2.2,2.7,3.2,3.7 \
    --rpts=6,6,6,10,15 \
    --bolus=1.5 \
    --slicedt=0.0185 \
    --fixbolus \
    --sliceband=6 \
    --fslanat=<anat_dir> \
    --pvcorr \
    -c <calibration_images> \
    --cmethod=voxel \
    --te=19 \
    --mc \
    --spatial=off \
    --fmap=<fmap> \
    --fmapmag=<averaged_corrected_fieldmap> \
    --fmapmagbrain=<averaged_corrected_fieldmap_brain> \
    --echospacing=0.00057 \
    --pedir=y \
    --gdcwarp=<gdc_warp>

### STEP 5
# Transform outputted perfusion results to native T1 surface space,
# and reference surface space (CIFTI)

