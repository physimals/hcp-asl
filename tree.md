```python
SUBJECTID_V1_MR
├── ASL #Contains intermediate pipeline outputs
│   ├── Calib #Calibration image: bias, distortion and banding correction
│   │   ├── Calib0 #PE calibration image
│   │   │   ├── BiasCorr
│   │   │   │   ├── SEbased
│   │   │   │   ├── calib0_bias.nii.gz
│   │   │   │   └── calib0_gdc_dc_restore.nii.gz #Bias and distortion corrected PE calibration image
│   │   │   ├── DistCorr
│   │   │   │   ├── asl2struct.mat #BBR Registration between PE calibration image and T1
│   │   │   │   ├── gdc_dc_calib0.nii.gz #Distortion corrected PE calibration image
│   │   │   │   └── struct2asl.mat
│   │   │   ├── MTCorr
│   │   │   │   ├── calib0_mtcorr_gdc_dc_restore.nii.gz #MT, distortion and bias corrected PE calibration image
│   │   │   │   └── mtcorr_gdc_calib0.nii.gz
│   │   │   ├── aslfs_mask.nii.gz
│   │   │   ├── calib0.nii.gz #Raw PE calibration image
│   │   │   └── calib_timing.nii.gz #Effective TI volume for oPE calibration image
│   │   ├── Calib1 #oPE calibration image#
│   │   │   ├── BiasCorr
│   │   │   │   ├── SEbased
│   │   │   │   ├── calib1_bias.nii.gz
│   │   │   │   └── calib1_gdc_dc_restore.nii.gz #Bias and distortion corrected oPE calibration image
│   │   │   ├── DistCorr
│   │   │   │   ├── asl2struct.mat  #BBR Registration between oPE calibration image and T1
│   │   │   │   ├── gdc_dc_calib1.nii.gz #Distortion corrected oPE calibration image
│   │   │   │   └── struct2asl.mat
│   │   │   ├── MTCorr
│   │   │   │   ├── calib1_mtcorr_gdc_dc_restore.nii.gz #MT, distortion and bias corrected oPE calibration image
│   │   │   │   └── mtcorr_gdc_calib1.nii.gz 
│   │   │   ├── aslfs_mask.nii.gz
│   │   │   └── calib1.nii.gz #Raw oPE calibration image
│   │   └── correct_M0.log
│   ├── TIs #ASL timeseries
│   │   ├── Betas
│   │   ├── BiasCorr
│   │   │   ├── tis_biascorr.nii.gz #Bias corrected ASL timeseries
│   │   │   └── tis_gdc_dc_restore.nii.gz #Bias and distortion corrected ASL timeseries
│   │   ├── DistCorr
│   │   │   └── FirstPass
│   │   │       └── tis_gdc.nii.gz #Gradient distortion corrected ASL timeseries
│   │   ├── MTCorr
│   │   │   ├── tis_gdc_dc_restore_mtcorr.nii.gz #MT and distortion corrected ASL timeseries
│   │   ├── MoCo #Motion correction
│   │   │   ├── asln2m0_final.mat #Motion correction matrices referenced to PE calibration image
│   │   │   ├── final_registration_TIs #ASL timeseries fully corrected and registered to PE calibration image
│   │   │   └── tis_gdc_dc_moco.nii.gz #Motion and distortion corrected ASL timeseries
│   │   ├── OxfordASL #BASIC perfusion estimates, used for registration
│   │   │   ├── native_space #Perfusion parameters, native space#
│   │   │   │   ├── aCBV.nii.gz
│   │   │   │   ├── arrival.nii.gz
│   │   │   │   ├── arrival_var.nii.gz
│   │   │   │   ├── mask.nii.gz
│   │   │   │   ├── perfusion.nii.gz
│   │   │   │   └── perfusion_var.nii.gz
│   │   │   └── oxford_asl.log
│   │   ├── STCorr
│   │   ├── STCorr2 #ST correction, final run
│   │   ├── SatRecov
│   │   ├── SatRecov2 #Saturation recovery fitting, final run
│   │   ├── aslfs_mask.nii.gz
│   │   ├── combined_scaling_factors.nii.gz #Combined banding correction factors
│   │   ├── tis.nii.gz #Raw ASL timeseries#
│   │   └── tis_gdc_dc_moco_restore_bandcorr.nii.gz #Fully corrected ASL timeseries in native acquisition space
│   ├── distortion_estimation.log
│   ├── gradient_unwarp #Gradient distortion correction
│   ├── single_step_resample_to_asl0.log
│   └── topup #Susecptibility distortion correction
├── MNINonLinear #Standard space results, follows structure of T1w directory
│   ├── ASL #Perfusion parameters in standard space#
│       ├── Results
│       │   ├── Native
│       │   │   └── T1wOutputtoCIFTI
│       │   │       ├── OutputtoCIFTI
│       │   │       │   ├── arrival.L.badvert_ribbonroi.native.func.gii
│       │   │       │   ├── arrival.L.native.func.gii
│       │   │       │   ├── arrival.R.badvert_ribbonroi.native.func.gii
│       │   │       │   ├── arrival.R.native.func.gii
│       │   │       │   └── arrival_precision.nii.gz
│       │   │       ├── perfusion_calib.L.badvert_ribbonroi.native.func.gii
│       │   │       ├── perfusion_calib.L.native.func.gii
│       │   │       ├── perfusion_calib.R.badvert_ribbonroi.native.func.gii
│       │   │       ├── perfusion_calib.R.native.func.gii
│       │   │       └── perfusion_calib_precision.nii.gz
│       │   ├── OutputtoCIFTI
│       │   │   ├── arrival_Atlas.dscalar.nii
│       │   │   ├── arrival_MNI.nii.gz
│       │   │   ├── asl_grid_mni.nii.gz
│       │   │   ├── perfusion_calib_Atlas.dscalar.nii
│       │   │   └── perfusion_calib_MNI.nii.gz
│       │   ├── OxfordASL
│       │   │   └── std_space
│       │   │       ├── arrival.nii.gz
│       │   │       ├── arrival_var.nii.gz
│       │   │       ├── perfusion_calib.nii.gz
│       │   │       ├── perfusion_var_calib.nii.gz
│       │   │       └── pvcorr
│       │   │           ├── arrival_gm_masked.nii.gz
│       │   │           ├── arrival_gm_var_masked.nii.gz
│       │   │           ├── arrival_wm_masked.nii.gz
│       │   │           ├── arrival_wm_var_masked.nii.gz
│       │   │           ├── perfusion_gm_calib_masked.nii.gz
│       │   │           ├── perfusion_gm_var_calib_masked.nii.gz
│       │   │           ├── perfusion_wm_calib_masked.nii.gz
│       │   │           └── perfusion_wm_var_calib_masked.nii.gz
│       │   └── pvcorr
│       │       ├── Native
│       │       │   └── T1wOutputtoCIFTI
│       │       │       ├── OutputtoCIFTI
│       │       │       │   ├── arrival.L.badvert_ribbonroi.native.func.gii
│       │       │       │   ├── arrival.L.native.func.gii
│       │       │       │   ├── arrival.R.badvert_ribbonroi.native.func.gii
│       │       │       │   ├── arrival.R.native.func.gii
│       │       │       │   └── arrival_precision.nii.gz
│       │       │       ├── perfusion_calib.L.badvert_ribbonroi.native.func.gii
│       │       │       ├── perfusion_calib.L.native.func.gii
│       │       │       ├── perfusion_calib.R.badvert_ribbonroi.native.func.gii
│       │       │       ├── perfusion_calib.R.native.func.gii
│       │       │       └── perfusion_calib_precision.nii.gz
│       │       └── OutputtoCIFTI
│       │           ├── arrival_Atlas.dscalar.nii
│       │           ├── arrival_MNI.nii.gz
│       │           ├── asl_grid_mni.nii.gz
│       │           ├── perfusion_calib_Atlas.dscalar.nii
│       │           └── perfusion_calib_MNI.nii.gz
│       ├── arrival_Atlas.dscalar.nii
│       ├── arrival_cifti_mean_nonzero.txt
│       ├── perfusion_calib_Atlas.dscalar.nii
│       ├── perfusion_calib_cifti_mean_nonzero.txt
│       ├── pvcorr_arrival_Atlas.dscalar.nii
│       ├── pvcorr_arrival_cifti_mean_nonzero.txt
│       ├── pvcorr_perfusion_calib_Atlas.dscalar.nii
│       └── pvcorr_perfusion_calib_cifti_mean_nonzero.txt
├── T1w #Subject anatomical results
│   ├── ACPCAlignment
│   ├── ASL #ASL-gridded T1w algined space (used for all perfusion analysis)
│   │   ├── Calib
│   │   │   └── Calib0
│   │   │       ├── DistCorr
│   │   │       ├── SEbased
│   │   │       ├── calib0_corr_aslt1w.nii.gz
│   │   │       ├── calib_aslt1w_stcorr_factors.nii.gz
│   │   │       └── calib_aslt1w_timing.nii.gz
│   │   ├── PVEs #Partial volume estimates
│   │   ├── Results #Surface projected perfusion for constructing CIFTIs
│   │   │   ├── OutputtoCIFTI #Without PVEc#
│   │   │   │   ├── arrival.L.badvert_ribbonroi.native.func.gii
│   │   │   │   ├── arrival.L.native.func.gii
│   │   │   │   ├── arrival.R.badvert_ribbonroi.native.func.gii
│   │   │   │   ├── arrival.R.native.func.gii
│   │   │   │   ├── arrival_precision.nii.gz
│   │   │   │   ├── perfusion_calib.L.badvert_ribbonroi.native.func.gii
│   │   │   │   ├── perfusion_calib.L.native.func.gii
│   │   │   │   ├── perfusion_calib.R.badvert_ribbonroi.native.func.gii
│   │   │   │   ├── perfusion_calib.R.native.func.gii
│   │   │   │   └── perfusion_calib_precision.nii.gz
│   │   │   └── pvcorr
│   │   │       └── OutputtoCIFTI #With PVEc
│   │   │           ├── arrival.L.badvert_ribbonroi.native.func.gii
│   │   │           ├── arrival.L.native.func.gii
│   │   │           ├── arrival.R.badvert_ribbonroi.native.func.gii
│   │   │           ├── arrival.R.native.func.gii
│   │   │           ├── arrival_precision.nii.gz
│   │   │           ├── perfusion_calib.L.badvert_ribbonroi.native.func.gii
│   │   │           ├── perfusion_calib.L.native.func.gii
│   │   │           ├── perfusion_calib.R.badvert_ribbonroi.native.func.gii
│   │   │           ├── perfusion_calib.R.native.func.gii
│   │   │           └── perfusion_calib_precision.nii.gz
│   │   ├── TIs
│   │   │   ├── Betas
│   │   │   ├── DistCorr
│   │   │   │   └── tis_gdc_dc_moco.nii.gz
│   │   │   ├── OxfordASL #Final perfusion estimation
│   │   │   │   ├── calib
│   │   │   │   │   ├── M0.txt
│   │   │   │   │   ├── logfile
│   │   │   │   │   └── refmask.nii.gz
│   │   │   │   ├── logfile
│   │   │   │   ├── native_space #Perfusion estimates
│   │   │   │   │   ├── aCBV.nii.gz
│   │   │   │   │   ├── aCBV_calib.nii.gz
│   │   │   │   │   ├── arrival.nii.gz
│   │   │   │   │   ├── arrival_gm_mean.txt
│   │   │   │   │   ├── arrival_var.nii.gz
│   │   │   │   │   ├── arrival_wm_mean.txt
│   │   │   │   │   ├── gm_mask.nii.gz
│   │   │   │   │   ├── gm_roi.nii.gz
│   │   │   │   │   ├── mask.nii.gz
│   │   │   │   │   ├── perfusion.nii.gz
│   │   │   │   │   ├── perfusion_calib.nii.gz
│   │   │   │   │   ├── perfusion_calib_gm_mean.txt
│   │   │   │   │   ├── perfusion_calib_wm_mean.txt
│   │   │   │   │   ├── perfusion_gm_mean.txt
│   │   │   │   │   ├── perfusion_norm.nii.gz
│   │   │   │   │   ├── perfusion_var.nii.gz
│   │   │   │   │   ├── perfusion_var_calib.nii.gz
│   │   │   │   │   ├── perfusion_var_norm.nii.gz
│   │   │   │   │   ├── perfusion_wm_mean.txt
│   │   │   │   │   ├── pvcorr #PVEc perfusion estimates#
│   │   │   │   │   │   ├── aCBV.nii.gz
│   │   │   │   │   │   ├── aCBV_calib.nii.gz
│   │   │   │   │   │   ├── arrival.nii.gz
│   │   │   │   │   │   ├── arrival_gm_mean.txt
│   │   │   │   │   │   ├── arrival_masked.nii.gz
│   │   │   │   │   │   ├── arrival_var.nii.gz
│   │   │   │   │   │   ├── arrival_var_masked.nii.gz
│   │   │   │   │   │   ├── arrival_wm.nii.gz
│   │   │   │   │   │   ├── arrival_wm_masked.nii.gz
│   │   │   │   │   │   ├── arrival_wm_var.nii.gz
│   │   │   │   │   │   ├── arrival_wm_var_masked.nii.gz
│   │   │   │   │   │   ├── arrival_wm_wm_mean.txt
│   │   │   │   │   │   ├── perfusion.nii.gz
│   │   │   │   │   │   ├── perfusion_calib.nii.gz
│   │   │   │   │   │   ├── perfusion_calib_gm_mean.txt
│   │   │   │   │   │   ├── perfusion_calib_masked.nii.gz
│   │   │   │   │   │   ├── perfusion_gm_mean.txt
│   │   │   │   │   │   ├── perfusion_masked.nii.gz
│   │   │   │   │   │   ├── perfusion_norm.nii.gz
│   │   │   │   │   │   ├── perfusion_var.nii.gz
│   │   │   │   │   │   ├── perfusion_var_calib.nii.gz
│   │   │   │   │   │   ├── perfusion_var_calib_masked.nii.gz
│   │   │   │   │   │   ├── perfusion_var_norm.nii.gz
│   │   │   │   │   │   ├── perfusion_wm.nii.gz
│   │   │   │   │   │   ├── perfusion_wm_calib.nii.gz
│   │   │   │   │   │   ├── perfusion_wm_calib_masked.nii.gz
│   │   │   │   │   │   ├── perfusion_wm_calib_wm_mean.txt
│   │   │   │   │   │   ├── perfusion_wm_masked.nii.gz
│   │   │   │   │   │   ├── perfusion_wm_norm.nii.gz
│   │   │   │   │   │   ├── perfusion_wm_var.nii.gz
│   │   │   │   │   │   ├── perfusion_wm_var_calib.nii.gz
│   │   │   │   │   │   ├── perfusion_wm_var_calib_masked.nii.gz
│   │   │   │   │   │   ├── perfusion_wm_var_norm.nii.gz
│   │   │   │   │   │   └── perfusion_wm_wm_mean.txt
│   │   │   │   │   ├── pvgm_inasl.nii.gz
│   │   │   │   │   ├── pvwm_inasl.nii.gz
│   │   │   │   │   ├── wm_mask.nii.gz
│   │   │   │   │   └── wm_roi.nii.gz
│   │   │   │   └── oxford_aslt1w.log
│   │   │   ├── asl_corr.nii.gz
│   │   │   ├── combined_scaling_factors.nii.gz
│   │   │   ├── reg #Registration between native ASL space and T1w aligned space
│   │   │   │   ├── ASL_FoV_brain_mask.nii.gz
│   │   │   │   ├── ASL_FoV_mask.nii.gz
│   │   │   │   ├── ASL_grid_T1w_brain_mask.nii.gz
│   │   │   │   ├── asl2orig_mgz_initial_bbr.dat
│   │   │   │   ├── asl2orig_mgz_initial_bbr.dat.log
│   │   │   │   ├── asl2orig_mgz_initial_bbr.dat.log.old
│   │   │   │   ├── asl2orig_mgz_initial_bbr.dat.mincost
│   │   │   │   ├── asl2orig_mgz_initial_bbr.dat.param
│   │   │   │   ├── asl2orig_mgz_initial_bbr.dat.sum
│   │   │   │   ├── asl2orig_mgz_initial_bbr.dat~
│   │   │   │   ├── asl2struct.mat
│   │   │   │   ├── asl2struct_reg.log
│   │   │   │   ├── fmapmag_aslt1w.nii.gz
│   │   │   │   └── mean_T1t_filt_aslt1w.nii.gz
│   │   │   ├── single_step_resample_to_aslt1w.log
│   │   │   └── timing_img_aslt1w.nii.gz
│   │   ├── aCBV_calib.nii.gz #These files are copies of those found in OxfordASL/native_space
│   │   ├── arrival.nii.gz
│   │   ├── arrival_gm_mean.txt
│   │   ├── arrival_var.nii.gz
│   │   ├── arrival_wm_mean.txt
│   │   ├── perfusion_calib.nii.gz
│   │   ├── perfusion_calib_gm_mean.txt
│   │   ├── perfusion_calib_wm_mean.txt
│   │   ├── perfusion_var_calib.nii.gz
│   │   ├── pvcorr_aCBV_calib.nii.gz
│   │   ├── pvcorr_arrival_gm_masked.nii.gz
│   │   ├── pvcorr_arrival_gm_mean.txt
│   │   ├── pvcorr_arrival_gm_var_masked.nii.gz
│   │   ├── pvcorr_arrival_wm_masked.nii.gz
│   │   ├── pvcorr_arrival_wm_mean.txt
│   │   ├── pvcorr_arrival_wm_var_masked.nii.gz
│   │   ├── pvcorr_perfusion_calib_gm_mean.txt
│   │   ├── pvcorr_perfusion_calib_wm_mean.txt
│   │   ├── pvcorr_perfusion_gm_calib_masked.nii.gz
│   │   ├── pvcorr_perfusion_gm_var_calib_masked.nii.gz
│   │   ├── pvcorr_perfusion_wm_calib_masked.nii.gz
│   │   ├── pvcorr_perfusion_wm_var_calib_masked.nii.gz
│   │   └── reg
│   │       ├── ASL_grid_T1w_acpc_dc_restore.nii.gz #Fully corrected T1 image in ASL-gridded space
│   │       └── wmmask.nii.gz
│   ├── BiasField_acpc_dc.nii.gz
│   ├── HCA7715581_V1_MR #FreeSurfer directory; surfaces and segmentations are used
│   ├── Native #Native resolution cortical surfaces onto which perfusion is projected
└── unprocessed #Contains raw pipeline inputs
    └── mbPCASLhr
        ├── HCA7715581_V1_MR_PCASLhr_SpinEchoFieldMap_AP.json
        ├── HCA7715581_V1_MR_PCASLhr_SpinEchoFieldMap_AP.nii.gz
        ├── HCA7715581_V1_MR_PCASLhr_SpinEchoFieldMap_PA.json
        ├── HCA7715581_V1_MR_PCASLhr_SpinEchoFieldMap_PA.nii.gz
        ├── HCA7715581_V1_MR_mbPCASLhr_PA.json
        ├── HCA7715581_V1_MR_mbPCASLhr_PA.nii.gz
```