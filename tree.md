```python
SUBJECTID_V1_MR
ASL # Native acquisition space processing 
├── Calib # Intermediate calibration image processing 
│   ├── Calib0
│   │   ├── BiasCorr
│   │   │   ├── SEbased
│   │   │   │   ├── AllGreyMatter.nii.gz
│   │   │   │   ├── [ .... ]
│   │   │   │   └── sebased_bias_dil.nii.gz
│   │   │   ├── gdc_dc_calib0_bias.nii.gz
│   │   │   └── gdc_dc_calib0_restore.nii.gz
│   │   ├── DistCorr
│   │   │   ├── asl2orig_mgz_initial_bbr.dat
│   │   │   ├── [ .... ]
│   │   │   └── struct2asl.mat
│   │   ├── MTCorr
│   │   │   ├── mtcorr_gdc_calib0.nii.gz
│   │   │   └── mtcorr_gdc_dc_calib0_restore.nii.gz
│   │   ├── aslfs_mask.nii.gz
│   │   ├── calib0.nii.gz
│   │   └── calib_timing.nii.gz
│   └── Calib1
│       ├── BiasCorr
│       │   ├── SEbased
│       │   │   ├── AllGreyMatter.nii.gz
│       │   │   ├── [ .... ]
│       │   │   └── sebased_bias_dil.nii.gz
│       │   ├── gdc_dc_calib1_bias.nii.gz
│       │   └── gdc_dc_calib1_restore.nii.gz
│       ├── DistCorr
│       │   ├── asl2orig_mgz_initial_bbr.dat
│       │   ├── [ .... ]
│       │   └── struct2asl.mat
│       ├── MTCorr
│       │   ├── mtcorr_gdc_calib1.nii.gz
│       │   └── mtcorr_gdc_dc_calib1_restore.nii.gz
│       ├── aslfs_mask.nii.gz
│       └── calib1.nii.gz
├── OxfordASL # Basic perfusion estimates for BBR registration 
│   ├── logfile
│   ├── native_space
│   │   ├── aCBV.nii.gz
│   │   ├── arrival.nii.gz
│   │   ├── arrival_var.nii.gz
│   │   ├── mask.nii.gz
│   │   ├── perfusion.nii.gz
│   │   └── perfusion_var.nii.gz
│   └── oxford_asl_inputs
│       ├── beta_perf.nii.gz
│       └── brain_fov_mask.nii.gz
├── TIs # Intermediate ASL timeseries processing 
│   ├── BiasCorr
│   │   ├── tis_biascorr.nii.gz
│   │   └── tis_dc_restore.nii.gz
│   ├── DistCorr
│   │   └── FirstPass
│   │       ├── dc_tis.nii.gz
│   │       ├── temp_reg_dc_tis.nii.gz
│   │       └── tis_gdc.nii.gz
│   ├── MTCorr
│   │   ├── tis_dc_restore_mtcorr.nii.gz
│   │   ├── tis_mtcorr.nii.gz
│   │   ├── tis_mtcorr_even.nii.gz
│   │   └── tis_mtcorr_odd.nii.gz
│   ├── MoCo
│   │   ├── asln2m0.mat
│   │   │   ├── MAT_0000
│   │   │   ├── [ .... ]
│   │   │   └── MAT_0085
│   │   ├── asln2m0_final.mat
│   │   │   ├── MAT_0000
│   │   │   ├── [ .... ]
│   │   │   └── MAT_0085
│   │   ├── final_registration_TIs.nii.gz
│   │   ├── final_registration_TIs.nii.gz.par
│   │   ├── fov_mask.nii.gz # Final motion-FoV mask 
│   │   ├── fov_mask_initial.nii.gz
│   │   ├── temp_reg_dc_tis_mtcorr_odd.nii.gz
│   │   └── tis_dc_moco.nii.gz
│   ├── MotionSubtraction
│   │   ├── beta_baseline.nii.gz
│   │   ├── beta_perf.nii.gz
│   │   ├── combined_mask.nii.gz
│   │   └── difference_mask.nii.gz
│   ├── STCorr
│   │   ├── st_scaling_factors.nii.gz
│   │   └── tis_stcorr.nii.gz
│   ├── STCorr2
│   │   ├── combined_scaling_factors_asln.nii.gz
│   │   ├── st_scaling_factors.nii.gz
│   │   └── tis_dc_restore_mtcorr_stcorr.nii.gz
│   ├── SatRecov
│   │   ├── nospatial
│   │   │   ├── finalMVN.nii.gz
│   │   │   ├── logfile
│   │   │   ├── mean_A.nii.gz
│   │   │   ├── mean_M0t.nii.gz
│   │   │   └── mean_T1t.nii.gz
│   │   └── spatial
│   │       ├── logfile
│   │       ├── mean_A.nii.gz
│   │       ├── mean_M0t.nii.gz
│   │       ├── mean_T1t.nii.gz
│   │       └── mean_T1t_filt.nii.gz
│   ├── SatRecov2
│   │   ├── nospatial
│   │   │   ├── finalMVN.nii.gz
│   │   │   ├── logfile
│   │   │   ├── mean_A.nii.gz
│   │   │   ├── mean_M0t.nii.gz
│   │   │   └── mean_T1t.nii.gz
│   │   └── spatial
│   │       ├── logfile
│   │       ├── mean_A.nii.gz
│   │       ├── mean_M0t.nii.gz
│   │       ├── mean_T1t.nii.gz
│   │       ├── mean_T1t_filt.nii.gz
│   │       └── mean_T1t_filt_asln.nii.gz
│   ├── brain_fov_mask.nii.gz # Combined brain and motion-FoV mask 
│   ├── brain_fov_mask_initial.nii.gz
│   ├── brain_mask.nii.gz
│   ├── combined_scaling_factors.nii.gz
│   ├── tis.nii.gz # Raw ASL timeseries 
│   ├── tis_dc_moco_restore.nii.gz
│   └── tis_dc_moco_restore_bandcorr.nii.gz # Fully corrected timeseries
├── gradient_unwarp
│   ├── fullWarp_abs.nii.gz
│   └── gdc_corr_vol1.nii.gz
├── topup
│   ├── Jacobian_01.nii.gz
│   ├── [ .... ]
│   ├── fmap_struct_reg
│   │   ├── asl2orig_mgz_initial_bbr.dat
│   │   ├── [ .... ]
│   │   └── asl2struct.mat
│   └── topup_params.txt
MNINonLinear # MNI space results 
├── ASL
│   ├── CIFTIPrepare
│   │   ├── arrival_Atlas.dscalar.nii
│   │   ├── [ .... ]
│   │   ├── perfusion_var_calib_MNI.nii.gz
│   │   └── pvcorr
│   │       ├── arrival_Atlas.dscalar.nii
│   │       ├── [ .... ]
│   │       └── perfusion_var_calib_MNI.nii.gz
│   ├── OxfordASL # Oxford ASL outputs (in MNI)
│   │   └── std_space
│   │       ├── arrival.nii.gz
│   │       ├── arrival_var.nii.gz
│   │       ├── perfusion_calib.nii.gz
│   │       ├── perfusion_var_calib.nii.gz
│   │       └── pvcorr
│   │           ├── arrival_gm_masked.nii.gz
│   │           ├── arrival_gm_var_masked.nii.gz
│   │           ├── arrival_wm_masked.nii.gz
│   │           ├── arrival_wm_var_masked.nii.gz
│   │           ├── perfusion_gm_calib_masked.nii.gz
│   │           ├── perfusion_gm_var_calib_masked.nii.gz
│   │           ├── perfusion_wm_calib_masked.nii.gz
│   │           └── perfusion_wm_var_calib_masked.nii.gz
│   ├── arrival.nii.gz # MNI NIFTI, ATT, non-PVEc
│   ├── arrival_Atlas.dscalar.nii # MNI CIFTI, ATT, non-PVEc
│   ├── arrival_cifti_mean_nonzero.txt # ROI calculation 
│   ├── perfusion_calib.nii.gz # MNI NIFTI, perfusion, non-PVEc
│   ├── perfusion_calib_Atlas.dscalar.nii # MNI CIFTI, perfusion, non-PVEc
│   ├── perfusion_calib_cifti_mean_nonzero.txt # ROI calculation 
│   ├── pvcorr_arrival_Atlas.dscalar.nii # MNI CIFTI, ATT, PVEc
│   ├── pvcorr_arrival_cifti_mean_nonzero.txt # ROI calculation 
│   ├── pvcorr_arrival_gm_masked.nii.gz # MNI NIFTI, ATT, PVEc GM
│   ├── pvcorr_arrival_wm_masked.nii.gz # MNI NIFTI, ATT, PVEc WM
│   ├── pvcorr_perfusion_calib_Atlas.dscalar.nii # MNI CIFTI, perfusion, PVEc
│   ├── pvcorr_perfusion_calib_cifti_mean_nonzero.txt # ROI calculation 
│   ├── pvcorr_perfusion_gm_calib_masked.nii.gz # MNI NIFTI, perfusion, PVEc GM
│   └── pvcorr_perfusion_wm_calib_masked.nii.gz # MNI NIFTI, perfusion, PVEc WM
T1w # T1w aligned space
└── ASL # T1w aligned at ASL resolution, used for final processing of all data. 
    ├── CIFTIPrepare
    │   ├── arrival.L.badvert_ribbonroi.native.func.gii
    │   ├── [ .... ]
    │   ├── perfusion_var_calib_s2.atlasroi.R.32k_fs_LR.func.gii
    │   └── pvcorr
    │       ├── arrival.L.badvert_ribbonroi.native.func.gii
    │       ├── [ .... ]
    │       └── perfusion_var_calib_s2.atlasroi.R.32k_fs_LR.func.gii
    ├── Calib
    │   └── Calib0
    │       ├── DistCorr
    │       │   ├── calib0_dc.nii.gz
    │       │   └── calib_mt_scaling_factors.nii.gz
    │       ├── SEbased
    │       │   ├── AllGreyMatter.nii.gz
    │       │   ├── [ .... ]
    │       │   └── sebased_bias_dilall.nii.gz
    │       ├── calib0_corr.nii.gz
    │       ├── calib_aslt1w_stcorr_factors.nii.gz
    │       └── calib_aslt1w_timing.nii.gz
    ├── HCD0378150_V1_MR_hcp_asl.log # HCP-ASL pipeline logfile 
    ├── OxfordASL
    │   ├── calib
    │   │   ├── M0.txt
    │   │   ├── logfile
    │   │   └── refmask.nii.gz
    │   ├── logfile
    │   ├── native_space
    │   │   ├── aCBV.nii.gz
    │   │   ├── [ .... ]
    │   │   ├── perfusion_wm_mean.txt
    │   │   ├── pvcorr
    │   │   │   ├── aCBV.nii.gz
    │   │   │   ├── [ .... ]
    │   │   │   └── perfusion_wm_wm_mean.txt
    │   │   └── wm_roi.nii.gz
    │   └── oxford_asl_inputs
    │       ├── beta_perf.nii.gz
    │       ├── brain_fov_mask.nii.gz
    │       ├── calib0_corr.nii.gz
    │       ├── pve_GM.nii.gz
    │       ├── pve_WM.nii.gz
    │       ├── timing_img_aslt1w.nii.gz
    │       └── vent_csf_mask.nii.gz
    ├── PVEs
    │   ├── pve_GM.nii.gz
    │   ├── pve_WM.nii.gz
    │   └── vent_csf_mask.nii.gz
    ├── TIs
    │   ├── DistCorr
    │   │   └── tis_dc_moco.nii.gz
    │   ├── MotionSubtraction
    │   │   ├── beta_baseline.nii.gz
    │   │   ├── beta_perf.nii.gz
    │   │   ├── combined_mask.nii.gz
    │   │   └── difference_mask.nii.gz
    │   ├── asl_corr.nii.gz
    │   ├── asl_noncorr.nii.gz
    │   ├── combined_scaling_factors.nii.gz
    │   ├── reg
    │   │   ├── ASL_grid_T1w_brain_mask.nii.gz
    │   │   ├── [ .... ]
    │   │   ├── asl2struct.mat
    │   │   ├── brain_fov_mask.nii.gz
    │   │   ├── fmapmag_aslt1w.nii.gz
    │   │   ├── fov_mask.nii.gz
    │   │   └── mean_T1t_filt_aslt1w.nii.gz
    │   └── timing_img_aslt1w.nii.gz
    ├── aCBV_calib.nii.gz # T1w NIFTI, arterial cerebral blood volume, non-PVEc
    ├── arrival.nii.gz # T1w NIFTI, ATT, non-PVEc
    ├── arrival_gm_mean.txt # ROI calculation 
    ├── arrival_var.nii.gz
    ├── arrival_wm_mean.txt # ROI calculation 
    ├── asl_corr.nii.gz # Fully corrected ASL timeseries 
    ├── asl_corr_subtracted.nii.gz # Fully corrected and subtracted ASL timeseries 
    ├── calib_corr.nii.gz # Fully corrected calibration image 
    ├── hcp_asl_report.ipynb # Jupyter notebook for QC 
    ├── perfusion_calib.nii.gz
    ├── perfusion_calib_gm_mean.txt
    ├── perfusion_calib_wm_mean.txt
    ├── perfusion_var_calib.nii.gz
    ├── pvcorr_aCBV_calib.nii.gz # T1w NIFTI, arterial cerebral blood volume, PVEc
    ├── pvcorr_arrival_gm_masked.nii.gz # T1w NIFTI, ATT, PVEc GM
    ├── pvcorr_arrival_gm_mean.txt # ROI calculation 
    ├── pvcorr_arrival_gm_var_masked.nii.gz
    ├── pvcorr_arrival_wm_masked.nii.gz # T1w NIFTI, ATT, PVEc WM
    ├── pvcorr_arrival_wm_mean.txt # ROI calculation 
    ├── pvcorr_arrival_wm_var_masked.nii.gz 
    ├── pvcorr_perfusion_calib_gm_mean.txt # ROI calculation 
    ├── pvcorr_perfusion_calib_wm_mean.txt # ROI calculation 
    ├── pvcorr_perfusion_gm_calib_masked.nii.gz # T1w NIFTI, perfusion, PVEc GM
    ├── pvcorr_perfusion_gm_var_calib_masked.nii.gz
    ├── pvcorr_perfusion_wm_calib_masked.nii.gz # T1w NIFTI, perfusion, PVEc WM
    ├── pvcorr_perfusion_wm_var_calib_masked.nii.gz
    ├── reg
    │   ├── ASL_grid_T1w_acpc_dc_restore.nii.gz
    │   └── wmmask.nii.gz
    └── roi_stats
        ├── fsl_identity.txt
        ├── roi_stats.csv
        ├── roi_stats_gm.csv
        └── roi_stats_wm.csv
```