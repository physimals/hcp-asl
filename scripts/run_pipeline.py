"""
This script performs the full minimal pre-processing ASL pipeline
for the Human Connectome Project (HCP) ASL data.

This currently requires that the script is called followed by
the directories of the subjects of interest and finally the
name of the empirical banding correction scaling factors image.
"""

import argparse
import logging
import os
from pathlib import Path
from shutil import copy, rmtree, copytree
from typing import Any
from dataclasses import dataclass

from hcpasl import __sha1__, __timestamp__, __version__
from hcpasl.asl_correction import initial_corrections_asl
from hcpasl.fully_corrected import fully_correct_asl_calibration_aslt1w
from hcpasl.asl_differencing import tag_control_differencing
from hcpasl.distortion_correction import derive_gdc_sdc
from hcpasl.key_outputs import copy_key_outputs
from hcpasl.calibration_correction import initial_corrections_calibration
from hcpasl.pv_estimation import run_pv_estimation
from hcpasl.qc import create_qc_report, roi_stats
from hcpasl.utils import (
    copy_oxford_asl_inputs,
    get_package_data_name,
    get_roi_stats_script,
    setup_logger,
    sp_run,
    split_asl,
    load_asl_params,
    AslParams,
    IBF,
    TIS,
    BOLUS,
    TE,
    SLICEDT,
    SLICEBAND,
    RPTS,
)


def process_subject(
    subid,
    subject_dir,
    eb_factors,
    mbpcasl,
    structural,
    fmaps,
    gradients,
    wmparc,
    ribbon,
    reg_name,
    territories_atlas,
    territories_labels,
    use_t1=False,
    cores=1,
    interpolation=3,
    nobandingcorr=False,
    outdir="hcp_asl",
    stages=set(range(14)),
    is_longitudinal=False,
    longitudinal_study_dir="",
    longitudinal_template="",
    asl_params: "AslParams" = None,
):
    """
    Run the hcp-asl pipeline for a given subject.

    Parameters
    ----------
    subid : str
        Subject ID.
    subject_dir : str
        Subject's data directory.
    eb_factors : pathlib.Path
        Path to a .txt file of pre-calculated empirical banding correction
        factors.
    mbpcasl : pathlib.Path
        Path to the subject's mbPCASL sequence.
    structural : dict
        Contains pathlib.Path locations of important structural
        files.
    fmaps : dict
        Contains pathlib.Path locations of the fieldmaps needed
        for distortion correction.
    gradients : str
        pathlib.Path to a gradient coefficients file for use in
        gradient distortion correction.
    wmparc : str
        pathlib.Path to wmparc.nii.gz from FreeSurfer for use in
        SE-based bias correction.
    ribbon : str
        pathlib.Path to ribbon.nii.gz from FreeSurfer for use in
        SE-based bias correction.
    reg_name : str
        Sphere to use for surface projection, e.g., MSMAll or MSMSulc.
    territories_atlas : pathlib.Path
        Path to vascular territories atlas.
    territories_labels: pathlib.Path
        Path to labels for a vascular territory atlas.
    use_t1 : bool, optional
        Whether or not to use the estimated T1 map in the
        oxford_asl run in structural space.
    cores : int, optional
        Number of cores to use.
        When applying motion correction, this is the number
        of cores that will be used by regtricks. Default is 1.
    interpolation : int, optional
        The interpolation order to use for registrations.
        Regtricks passes this on to scipy's map_coordinates.
        The meaning of the value can be found in the scipy
        documentation. Default is 3.
    nobandingcorr : bool, optional
        If this is True, the banding correction options in the
        pipeline will be switched off. Default is False (i.e.
        banding corrections are applied by default).
    outdir : str, optional
        Name of the results directory with subject folder.
    stages: set, optional
        List or set of integer stage numbers (zero-indexed) to run.
        All prior stages are assumed to have run successfully.
    is_longitudinal : bool, optional
        longitudinal mode
    longitudinal_study_dir: pathlib.Path
        Study dir that contains longitudinal sessions, only used in longitudinal mode
    longitudinal_template : str, required if is_longitudinal=True
        Longitudinal base template label, should match the one used
        in other HCP pipelines.

    """

    if not isinstance(stages, (list, set)):
        raise RuntimeError("stages must either be a set or list of ints")

    # copy the cross-sectional results to longitudinal folders
    @dataclass
    class RuntimeConfig:
        study_dir: Path = None
        subid: Any = None
        subject_dir: Path = None
        asl_dir: Path = None
        aslt1w_dir: Path = None
        label_control_dir: Path = None
        calib0_dir: Path = None
        calib1_dir: Path = None
        tis_name: Any = None
        calib0_name: Any = None
        calib1_name: Any = None
        gradunwarp_dir: Path = None
        topup_dir: Path = None
        gd_corr: Any = None
        t1w_dir: Path = None
        wmparc: Any = None
        ribbon: Any = None
        calib_corr: Any = None
        bias_field: Any = None
        calib2struct: Any = None
        asl_lc: Any = None
        scaling_factors: Any = None
        asl_subtract: Any = None
        asl0_brainmask: Any = None
        oxford_asl_dir: Path = None
        asl_scaling_factors: Any = None
        mt_name: Any = None
        t1_est: Any = None
        series: Any = None
        scaling_factors: Any = None

    confCross = RuntimeConfig(
        subject_dir=Path(subject_dir),
        subid=subid,
        study_dir=os.path.dirname(subject_dir),
    )
    confLong = RuntimeConfig(study_dir=longitudinal_study_dir)

    def assign_vars0(conf):
        conf.asl_dir, conf.aslt1w_dir = [
            conf.subject_dir / outdir / name for name in ("ASL", "T1w/ASL")
        ]
        conf.label_control_dir, conf.calib0_dir, conf.calib1_dir = [
            conf.asl_dir / name
            for name in ("label_control", "calibration/calib0", "calibration/calib1")
        ]
        for d in [
            conf.asl_dir,
            conf.aslt1w_dir,
            conf.label_control_dir,
            conf.calib0_dir,
            conf.calib1_dir,
        ]:
            d.mkdir(exist_ok=True, parents=True)

        # split ASL sequence into label-control label_control and calibration images
        conf.tis_name, conf.calib0_name, conf.calib1_name = [
            d / name
            for d, name in zip(
                (conf.label_control_dir, conf.calib0_dir, conf.calib1_dir),
                ("label_control.nii.gz", "calib0.nii.gz", "calib1.nii.gz"),
            )
        ]

    if is_longitudinal:
        if longitudinal_template == "":
            logging.error(
                "Longitudinal template label may not be empty in longitudinal mode."
            )
        if {0, 1, 2, 3, 4, 5} & stages:
            logging.error("Stages 0-5 cannot be run in longitudinal mode")
        # if 6 not in stages:
        #    logging.error("Stage 6 is required in longitudinal mode")
        confLong.subid = f"{subid}.long.{longitudinal_template}"
        confLong.subject_dir = Path(f"{confLong.study_dir}/{confLong.subid}")
        conf = confLong
        assign_vars0(confLong)
    else:
        conf = confCross

    assign_vars0(confCross)

    if 0 in stages:
        logging.info(
            "Stage 0: Splitting ASL sequence into label-control pairs and calibration images."
        )
        split_asl(
            mbpcasl,
            conf.tis_name,
            conf.calib0_name,
            conf.calib1_name,
            params=asl_params,
        )

    # Run gradient_unwarp and topup on the calibration images
    # results stored in gradunwarp_dir and topup_dir respectively

    def assign_vars1(conf):
        conf.gradunwarp_dir = conf.asl_dir / "gradient_unwarp"
        conf.topup_dir = conf.asl_dir / "topup"
        conf.gd_corr = gradients is not None

    assign_vars1(confCross)
    if is_longitudinal:
        assign_vars1(confLong)

    if 1 in stages:
        logging.info(
            "Stage 1: Derive gradient and susceptibility distortion correction."
        )
        if not conf.gd_corr:
            logging.info(
                "Gradient coefficient file not provided, derivation will be skipped."
            )
            logging.info("Deriving susceptibility distortion correction only.")
        derive_gdc_sdc(
            vol=str(conf.calib0_name),
            coeffs_path=gradients,
            gradunwarp_dir=conf.gradunwarp_dir,
            topup_dir=conf.topup_dir,
            pa_sefm=str(fmaps["PA"]),
            ap_sefm=str(fmaps["AP"]),
            interpolation=interpolation,
            gd_corr=conf.gd_corr,
        )

    # Apply corrections derived thus far to M0 image
    # Estimate biasfield
    def assign_vars2(conf, is_cross):
        if is_cross:
            conf.t1w_dir = structural["struct"].parent
            conf.wmparc = wmparc
            conf.ribbon = ribbon
        else:
            conf.t1w_dir = Path(f"{conf.subject_dir}/T1w")
            conf.wmparc = Path(f"{conf.t1w_dir}/wmparc.nii.gz")
            conf.ribbon = Path(f"{conf.t1w_dir}/ribbon.nii.gz")

    assign_vars2(confCross, True)
    if is_longitudinal:
        assign_vars2(confLong, False)

    if 2 in stages:
        logging.info("Stage 2: Derive and apply initial corrections to M0 image.")
        initial_corrections_calibration(
            subject_id=conf.subid,
            calib_dir=conf.calib0_dir.parent,
            eb_factors=eb_factors,
            t1w_dir=conf.t1w_dir,
            aslt1w_dir=conf.aslt1w_dir,
            gradunwarp_dir=conf.gradunwarp_dir,
            gd_corr=conf.gd_corr,
            topup_dir=conf.topup_dir,
            wmparc=conf.wmparc,
            ribbon=conf.ribbon,
            interpolation=interpolation,
            nobandingcorr=nobandingcorr,
        )

    def assign_vars3(conf):
        conf.calib_corr = conf.calib0_dir / "calib0_initial_corrected.nii.gz"
        conf.bias_field = conf.calib0_dir / "bias_correction/calib0_biasfield.nii.gz"
        conf.calib2struct = conf.calib0_dir / "registration/asl2struct.mat"

    assign_vars3(confCross)
    if is_longitudinal:
        assign_vars3(confLong)

    # Apply corrections derived thus far to ASL timeseries
    if 3 in stages:
        logging.info("Stage 3: Derive and apply initial corrections to ASL timeseries.")
        initial_corrections_asl(
            subject_dir=conf.subject_dir,
            label_control_dir=conf.label_control_dir,
            eb_factors=eb_factors,
            bias_name=conf.bias_field,
            calib_name=conf.calib_corr,
            calib2struct=conf.calib2struct,
            gradunwarp_dir=conf.gradunwarp_dir,
            gd_corr=conf.gd_corr,
            topup_dir=conf.topup_dir,
            t1w_dir=conf.t1w_dir,
            cores=cores,
            interpolation=interpolation,
            nobandingcorr=nobandingcorr,
            params=asl_params,
        )

    # perform tag-control subtraction in ASL0 space
    def assign_vars4(conf):
        conf.asl_lc = conf.label_control_dir / "label_control_corrected.nii.gz"
        conf.scaling_factors = (
            conf.label_control_dir / "combined_scaling_factors_mc.nii.gz"
        )
        conf.asl_subtract = conf.label_control_dir / "motion_subtraction"
        conf.asl0_brainmask = conf.label_control_dir / "brain_fov_mask.nii.gz"

    assign_vars4(confCross)
    if is_longitudinal:
        assign_vars4(confLong)

    if 4 in stages:
        logging.info("Stage 4: Label-control subtraction in native ASL space.")
        tag_control_differencing(
            conf.asl_lc,
            conf.scaling_factors,
            conf.asl_subtract,
            mask=conf.asl0_brainmask,
        )

    # estimate perfusion in ASL0 space using oxford_asl
    def assign_vars5(conf):
        conf.oxford_asl_dir = conf.label_control_dir.parent / "perfusion_estimation"

    assign_vars5(confCross)
    if is_longitudinal:
        assign_vars5(confLong)

    if 5 in stages:
        logging.info("Stage 5: Perfusion estimation in ASL native space.")
        logging.info(
            f"Copying oxford_asl inputs to one location ({str(conf.oxford_asl_dir / 'oxford_asl_inputs')})."
        )
        oxasl_inputs = {
            "-i": conf.asl_subtract / "beta_perf.nii.gz",
            "-m": conf.asl0_brainmask,
        }
        copy_oxford_asl_inputs(oxasl_inputs, conf.oxford_asl_dir / "oxford_asl_inputs")
        # use resolved parameters for the oxford_asl call
        oxford_asl_call = [
            "oxford_asl",
            *[f"{k}={str(v)}" for k, v in oxasl_inputs.items()],
            f"-o={str(conf.oxford_asl_dir)}",
            "--tis="
            + ",".join([str(t) for t in (asl_params.tis if asl_params else TIS)]),
            f"--slicedt={(asl_params.slicedt if asl_params else SLICEDT)}",
            f"--sliceband={(asl_params.sliceband if asl_params else SLICEBAND)}",
            "--casl",
            f"--ibf={(asl_params.ibf if asl_params else IBF)}",
            "--iaf=diff",
            "--rpts="
            + ",".join([str(r) for r in (asl_params.rpts if asl_params else RPTS)]),
            "--fixbolus",
            f"--bolus={(asl_params.bolus if asl_params else BOLUS)}",
            f"--te={(asl_params.te_ms if asl_params else TE)}",
            "--spatial=off",
            "--debug",
        ]
        if use_t1:
            est_t1 = (
                conf.label_control_dir
                / "saturation_recovery/second/spatial/mean_T1t_filt.nii.gz"
            )
            oxford_asl_call.append(f"--t1im={str(est_t1)}")
        logging.info(oxford_asl_call)
        sp_run(oxford_asl_call)

    # get data in ASLT1w space
    def assign_vars6(conf):
        if not nobandingcorr:
            conf.asl_scaling_factors = (
                conf.label_control_dir
                / "slicetime_correction/second/combined_scaling_factors_asln.nii.gz"
            )
            conf.mt_name = eb_factors
        else:
            conf.asl_scaling_factors, conf.mt_name = None, None
        conf.t1_est = (
            conf.label_control_dir
            / "saturation_recovery/second/spatial/mean_T1t_filt.nii.gz"
        )

    assign_vars6(confCross)
    if is_longitudinal:
        assign_vars6(confLong)
        copydirs = [
            "calibration",
            "gradient_unwarp",
            "label_control",
            "topup",
            "perfusion_estimation",
        ]
        # copy over the results of all previous stages.
        for dir in copydirs:
            if (confLong.asl_dir / dir).exists():
                rmtree(confLong.asl_dir / dir)
            subdir = confCross.asl_dir / dir
            if subdir.exists():
                copytree(subdir, confLong.asl_dir / dir, True)
        # this will need to be re-generated at stage 6.
        fmap_struct_reg = confLong.topup_dir / "fmap_struct_reg/asl2struct.mat"
        if fmap_struct_reg.exists():
            os.remove(fmap_struct_reg)

    # In longitudinal mode, this stage must run first. Stages 0-5 are skipped, with
    # results copied from corresponding cross-sectional folders.
    if 6 in stages:
        logging.info(
            "Stage 6: Fully-correct ASL and calibration into ASL-gridded T1w space."
        )
        fully_correct_asl_calibration_aslt1w(
            asl_name=conf.tis_name,
            calib_name=conf.calib0_name,
            subid=conf.subid,
            subject_dir=conf.subject_dir,
            t1w_dir=conf.t1w_dir,
            aslt1w_dir=conf.aslt1w_dir,
            moco_dir=conf.label_control_dir
            / "motion_correction/asln2calibration_final.mat",
            perfusion_name=conf.label_control_dir.parent
            / "perfusion_estimation/native_space/perfusion.nii.gz",
            gradunwarp_dir=conf.gradunwarp_dir,
            gd_corr=conf.gd_corr,
            topup_dir=conf.topup_dir,
            ribbon=conf.ribbon,
            wmparc=conf.wmparc,
            asl_scaling_factors=conf.asl_scaling_factors,
            eb_factors=conf.mt_name,
            t1_est=conf.t1_est,
            nobandingcorr=nobandingcorr,
            interpolation=interpolation,
            cores=cores,
            is_longitudinal=is_longitudinal,
            aslt1w_cross_dir=confCross.aslt1w_dir,
            topup_cross_dir=confCross.topup_dir,
        )
        copy(
            conf.aslt1w_dir / "calibration/calib0/calib0_corrected.nii.gz",
            conf.aslt1w_dir / "calib_corrected.nii.gz",
        )
        copy(
            conf.aslt1w_dir / "label_control/label_control_corrected.nii.gz",
            conf.aslt1w_dir / "label_control_corrected.nii.gz",
        )

    # starting with stage 7, all outputs in longitudinal mode are re-generated,
    # so no copying from cross-sectional to longitudinal occurs and multiple config structures aren't necessary.

    # perform partial volume estimation
    if 7 in stages:
        logging.info("Stage 7: Partial volume estimation in ASLT1w space.")
        run_pv_estimation(conf.subject_dir, cores, outdir, interpolation)

    # perform tag-control subtraction in ASLT1w space
    # aslt1w_dir = aslt1w_dir
    series = conf.aslt1w_dir / "label_control/label_control_corrected.nii.gz"
    scaling_factors = (
        conf.aslt1w_dir / "label_control/label_control_scaling_factors.nii.gz"
    )
    subtracted_dir = conf.aslt1w_dir / "label_control/motion_subtraction"
    brainmask = conf.aslt1w_dir / "registration/brain_fov_mask.nii.gz"

    if 8 in stages:
        logging.info("Stage 8: Label-control subtraction in ASLT1w space")
        tag_control_differencing(
            series, scaling_factors, subtracted_dir, mask=brainmask
        )
        copy(
            subtracted_dir / "beta_perf.nii.gz",
            conf.aslt1w_dir / "label_control_corrected_subtracted.nii.gz",
        )

    # final perfusion estimation in ASLT1w space
    pve_dir = conf.aslt1w_dir / "pvs"
    gm_pve, wm_pve = [pve_dir / f"pv_{tiss}.nii.gz" for tiss in ("GM", "WM")]
    oxford_aslt1w_dir = conf.aslt1w_dir / "perfusion_estimation"
    oxford_aslt1w_dir.mkdir(parents=True, exist_ok=True)
    if 9 in stages:
        logging.info("Stage 9: Perfusion estimation in ASLT1w space")
        logging.info(
            f"Copying oxford_asl inputs to one location ({str(oxford_aslt1w_dir / 'oxford_asl_inputs')})."
        )
        oxasl_inputs = {
            "-i": subtracted_dir / "beta_perf.nii.gz",
            "--pvgm": gm_pve,
            "--pvwm": wm_pve,
            "--csf": pve_dir / "vent_csf_mask.nii.gz",
            "-c": conf.aslt1w_dir / "calibration/calib0/calib0_corrected.nii.gz",
            "-m": conf.aslt1w_dir / "registration/brain_fov_mask.nii.gz",
            "--tiimg": conf.aslt1w_dir / "label_control/timing_img.nii.gz",
        }
        if use_t1:
            oxasl_inputs["--t1im"] = (
                conf.aslt1w_dir / "registration/mean_T1t_filt_aslt1w.nii.gz"
            )
        copy_oxford_asl_inputs(oxasl_inputs, oxford_aslt1w_dir / "oxford_asl_inputs")

        oxford_aslt1w_call = [
            "oxford_asl",
            f"-o={oxford_aslt1w_dir}",
            *[f"{k}={str(v)}" for k, v in oxasl_inputs.items()],
            "--casl",
            f"--ibf={IBF}",
            "--iaf=diff",
            "--rpts=" + ",".join([str(r) for r in RPTS]),
            "--fixbolus",
            f"--bolus={BOLUS}",
            f"--te={TE}",
            "--spatial=off",
            "--tr=8",
            "--pvcorr",
            "--debug",
        ]
        if use_t1:
            est_t1 = conf.aslt1w_dir / "registration/mean_T1t_filt.nii.gz"
            oxford_aslt1w_call.append(f"--t1im={str(est_t1)}")
        sp_run(oxford_aslt1w_call)

    mninonlinear_name = conf.subject_dir / "MNINonLinear"
    struct_name = (
        structural["struct"]
        if not is_longitudinal
        else Path(f"{conf.t1w_dir}/T1w_acpc_dc_restore.nii.gz")
    )

    if 10 in stages:
        logging.info("Stage 10: Summary statistics within ROIs.")
        roi_stats(
            struct_name=struct_name,
            oxford_asl_dir=oxford_aslt1w_dir,
            gm_pve=gm_pve,
            wm_pve=wm_pve,
            std2struct_name=mninonlinear_name / "xfms/standard2acpc_dc.nii.gz",
            roi_stats_dir=conf.aslt1w_dir / "roi_stats",
            territories_atlas=territories_atlas,
            territories_labels=territories_labels,
        )

    if 11 in stages:
        logging.info("Stage 11: Volume to surface projection.")
        surface_projection_stage(
            subject_dir=conf.subject_dir,
            subject_id=conf.subid,
            outdir=outdir,
            reg_name=reg_name,
        )

    if 12 in stages:
        logging.info(
            "Stage 12: Copy key results into $outdir/T1w/ASL and $outdir/MNINonLinear/ASL"
        )
        copy_outputs(conf.subject_dir, outdir)

    if 13 in stages:
        logging.info("Stage 13: Create QC workbench scene.")
        create_qc_report(
            subject_id=conf.subid,
            subject_dir=conf.subject_dir,
            outdir=outdir,
            reg_name=reg_name,
        )
    logging.info("Pipeline complete.")


def surface_projection_stage(
    subject_id,
    subject_dir,
    outdir,
    lowresmesh="32",
    FinalASLRes="2.5",
    SmoothingFWHM="2",
    GreyOrdsRes="2",
    reg_name="MSMAll",
):
    """
    Project perfusion results to the cortical surface and generate
    CIFTI representation which includes both low res mesh surfaces
    in MSMAll Atlas space, and subcortical structures in MNI
    voxel space

    Parameters
    ----------
    subject_dir : pathlib.Path
        Subject's data directory
    """

    # Projection scripts path:
    script = "PerfusionCIFTIProcessingPipelineASL.sh"
    wb_path = os.environ["CARET7DIR"]

    if not outdir:
        outdir = subject_dir
    else:
        outdir = subject_dir / outdir

    ASLVariable = ["perfusion_calib", "arrival", "perfusion_var_calib", "arrival_var"]
    ASLVariableVar = [
        "perfusion_var_calib",
        "arrival_var",
        "perfusion_var_calib",
        "arrival_var",
    ]

    for idx in range(4):
        non_pvcorr_cmd = [
            script,
            str(subject_dir),
            subject_id,
            ASLVariable[idx],
            ASLVariableVar[idx],
            lowresmesh,
            FinalASLRes,
            SmoothingFWHM,
            GreyOrdsRes,
            reg_name,
            wb_path,
            "false",
            str(outdir),
        ]

        pvcorr_cmd = [
            script,
            str(subject_dir),
            subject_id,
            ASLVariable[idx],
            ASLVariableVar[idx],
            lowresmesh,
            FinalASLRes,
            SmoothingFWHM,
            GreyOrdsRes,
            reg_name,
            wb_path,
            "true",
            str(outdir),
        ]

        sp_run(non_pvcorr_cmd)
        sp_run(pvcorr_cmd)


def copy_outputs(subject_dir, outdir):
    """
    Copy key pipeline outputs to the T1w and MNI aligned high level ASL directory

    Parameters
    ----------
    subject_dir : pathlib.Path
        Path to the subject data directory
    """

    path_to_outs = str(subject_dir / outdir)
    mni_raw = str(subject_dir / "MNINonLinear")
    t1w_preproc = str(subject_dir / "T1w")
    copy_key_outputs(path_to_outs, t1w_preproc, mni_raw)


def main():
    """
    Main entry point for the hcp-asl pipeline.
    """

    env_var = ["HCPPIPEDIR", "FREESURFER_HOME", "FSLDIR", "CARET7DIR"]
    for ev in env_var:
        if not bool(os.environ.get(ev)):
            raise RuntimeError(
                f"Environment variable {ev} must be set (see installation instructions)"
            )

    # Try and load the ROI stats script now - func will raise exception if not found.
    get_roi_stats_script()

    # argument handling
    parser = argparse.ArgumentParser(
        description="Minimal processing pipeline for HCP Lifespan ASL data."
    )

    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "--subid", help="Subject ID, including '_V1_MR'", required=True
    )
    required.add_argument(
        "--subdir",
        help="Subject's structural pre-processed data directory",
        required=True,
    )
    required.add_argument(
        "--mbpcasl",
        help="Filename for the mbPCASLhr acquisition",
        required=True,
    )
    required.add_argument(
        "--fmap_ap",
        help="Filename for the AP fieldmap for use in distortion correction",
        required=True,
    )
    required.add_argument(
        "--fmap_pa",
        help="Filename for the PA fieldmap for use in distortion correction",
        required=True,
    )
    optional = parser.add_argument_group(
        "optional arguments",
        description="(will attempt to load from default locations in $subdir)",
    )
    optional.add_argument(
        "--grads",
        help="Filename of the gradient coefficients for gradient"
        + " distortion correction.",
        required=False,
    )
    optional.add_argument(
        "--struct",
        help="Filename for the acpc-aligned, dc and restored structural image,"
        + " default is within subject's directory",
    )
    optional.add_argument(
        "--sbrain",
        help="Filename for the brain-extracted acpc-aligned, "
        + "dc-bcd structural image, default is within subject's directory",
    )
    optional.add_argument(
        "--wmparc",
        help="wmparc.nii.gz from FreeSurfer for use in SE-based bias correction",
    )
    optional.add_argument(
        "--ribbon",
        help="ribbon.nii.gz from FreeSurfer for use in SE-based bias correction,"
        + " default is within subject's directory",
    )
    optional.add_argument(
        "--regname",
        help="Sphere to use for surface projection stages, default is MSMAll (requires fMRI pre-processing), located at ${StudyFolder}/${Subject}/MNINonLinear/Native/${Subject}.${Hemisphere}.sphere.MSMSulc.native.surf.gii",
        default="MSMAll",
        dest="reg_name",
    )
    optional.add_argument(
        "--use_t1",
        help="If this flag is provided, the T1 estimates from the satrecov "
        + "will also be registered to ASL-gridded T1 space for use in "
        + "perfusion estimation via oxford_asl.",
        action="store_true",
    )
    optional.add_argument(
        "--mtname",
        help="Filename for scaling factors used for empirical banding  "
        + "correction. If not provided, the pipeline will "
        + "use the scaling factors included with the distribution.",
    )
    optional.add_argument(
        "--stages",
        help="Pipeline stages (zero-indexed, separated by spaces) to run, eg 0 3 5",
        nargs="+",
        type=int,
        default=None,
        metavar="N",
    )
    optional.add_argument(
        "--cores",
        help="Number of cores to use for multi-core operations. Default is 1.",
        default=1,
        type=int,
    )
    optional.add_argument(
        "--interpolation",
        help="Interpolation order for registrations. This can be any "
        + "integer from 0-5 inclusive. Default is 3. See scipy's "
        + "map_coordinates for more details.",
        default=3,
        type=int,
        choices=range(0, 5 + 1),
    )
    optional.add_argument(
        "--nobandingcorr",
        help="Don't apply empirical banding and slice-time corrections.",
        action="store_true",
    )
    optional.add_argument(
        "--territories_atlas",
        help="Location of vascular territory atlas.",
        default=get_package_data_name("vascular_territories_atlas.nii.gz"),
    )
    optional.add_argument(
        "--territories_labels",
        help="Location of txt file with labels for vascular territory atlas.",
        default=get_package_data_name("vascular_territories_atlas_labels.txt"),
    )
    optional.add_argument(
        "--outdir",
        help="Name of output directory to be created within subjects's directory.",
        default="",
    )
    optional.add_argument(
        "--clean",
        help="Remove all previous outputs found within --outdir",
        action="store_true",
    )
    optional.add_argument(
        "--is_longitudinal",
        help="Longitudinal processing",
        action="store_true",
    )
    optional.add_argument(
        "--longitudinal_study_dir",
        help="Study dir that contains longitudinal sessions, only used in longitudinal mode",
    )

    optional.add_argument(
        "--longitudinal_template",
        help="Longitudinal base template label (required in longitudinal mode)",
        default="",
    )

    # ASL acquisition parameter overrides
    optional.add_argument("--ntis", type=int, help="Number of TIs")
    optional.add_argument(
        "--tis",
        type=float,
        nargs="+",
        metavar="TI",
        help="Space-separated list of TIs in seconds (e.g., 1.7 2.2 2.7 3.2 3.7)",
    )
    optional.add_argument(
        "--rpts",
        type=int,
        nargs="+",
        metavar="RPT",
        help="Space-separated repeats for each TI (e.g., 6 6 6 10 15)",
    )
    optional.add_argument("--bolus", type=float, help="Labeling/bolus duration (s)")
    optional.add_argument("--slicedt", type=float, help="Slice time (s)")
    optional.add_argument(
        "--sliceband",
        type=int,
        help="Slices per band (if omitted, derived from sidecar MB factor)",
    )
    optional.add_argument(
        "--te",
        type=float,
        dest="te_ms",
        help="Echo time in milliseconds",
    )
    optional.add_argument(
        "--tail_discard_vols",
        type=int,
        help="Volumes immediately before calibrations to discard (default 2)",
    )
    optional.add_argument("--ibf", type=str, help="Input block format (default 'tis')")

    # assign arguments to variables
    args = parser.parse_args()
    subid = args.subid
    subject_dir = Path(args.subdir).resolve(strict=True)
    base_dir = subject_dir / Path(args.outdir)

    if args.clean:
        for d in ["ASL", "T1w/ASL", "MNINonLinear/ASL"]:
            if (base_dir / d).exists():
                rmtree(base_dir / d, ignore_errors=True)

    base_dir.mkdir(exist_ok=True, parents=True)
    log_path = base_dir / f"T1w/ASL/{subid}_hcp_asl.log"
    log_path.parent.mkdir(exist_ok=True, parents=True)
    setup_logger(log_path)

    logging.info(
        f"HCP-ASL pipeline v{__version__} (commit {__sha1__} on {__timestamp__})."
    )
    logging.info(f"Logging to {log_path}")

    # Look for required files in default paths if not provided.
    if args.struct is None:
        args.struct = subject_dir / "T1w/T1w_acpc_dc_restore.nii.gz"
        logging.info(f"Using default for struct: {args.struct}")
    if not os.path.exists(args.struct):
        raise ValueError(f"Path to struct does not exist: {args.struct}")

    if args.sbrain is None:
        args.sbrain = subject_dir / "T1w/T1w_acpc_dc_restore_brain.nii.gz"
        logging.info(f"Using default for sbrain: {args.sbrain}")
    if not os.path.exists(args.sbrain):
        raise ValueError(f"Path to sbrain does not exist: {args.sbrain}")

    if args.wmparc is None:
        args.wmparc = subject_dir / "T1w/wmparc.nii.gz"
        logging.info(f"Using default for wmparc: {args.wmparc}")
    if not os.path.exists(args.wmparc):
        raise ValueError(f"Path to wmparc does not exist: {args.wmparc}")

    if args.ribbon is None:
        args.ribbon = subject_dir / "T1w/ribbon.nii.gz"
        logging.info(f"Using default for ribbon: {args.ribbon}")
    if not os.path.exists(args.ribbon):
        raise ValueError(f"Path to ribbon does not exist: {args.ribbon}")

    # parse remaining arguments
    if args.mtname:
        mtname = Path(args.mtname).resolve(strict=True)
    elif not args.nobandingcorr:
        mtname = get_package_data_name("empirical_banding_factors.txt")
    else:
        mtname = None
    structural = {
        "struct": Path(args.struct).resolve(strict=True),
        "sbrain": Path(args.sbrain).resolve(strict=True),
    }
    mbpcasl = Path(args.mbpcasl).resolve(strict=True)

    # resolve ASL acquisition parameters
    asl_params: AslParams = load_asl_params(
        mbpcasl,
        ntis=args.ntis,
        tis=args.tis,
        rpts=args.rpts,
        bolus=args.bolus,
        slicedt=args.slicedt,
        te_ms=args.te_ms,
        sliceband=args.sliceband,
        tail_discard_vols=args.tail_discard_vols,
        ibf=args.ibf,
    )

    fmaps = {
        "AP": Path(args.fmap_ap).resolve(strict=True),
        "PA": Path(args.fmap_pa).resolve(strict=True),
    }
    if args.grads is not None:
        grads = Path(args.grads).resolve(strict=True)
    else:
        logging.info(
            "No gradient coefficients provided. Gradient distortion correction won't be performed."
        )
        grads = None

    if args.stages is None:
        if args.is_longitudinal:
            stages = set(range(6, 14))
        else:
            stages = set(range(14))
    else:
        stages = set(args.stages)

    longitudinal_study_dir = (
        None
        if args.longitudinal_study_dir is None
        else Path(args.longitudinal_study_dir)
    )

    logging.info("All pipeline arguments:")
    for k, v in vars(args).items():
        logging.info(f"{k}: {v}")
    logging.info(f"Resolved ASL params: {asl_params}")

    # process subject
    logging.info(f"Processing subject {subject_dir}.")
    try:
        process_subject(
            subid=subid,
            subject_dir=subject_dir,
            eb_factors=mtname,
            cores=args.cores,
            interpolation=args.interpolation,
            gradients=grads,
            mbpcasl=mbpcasl,
            territories_atlas=args.territories_atlas,
            territories_labels=args.territories_labels,
            structural=structural,
            fmaps=fmaps,
            use_t1=args.use_t1,
            wmparc=args.wmparc,
            ribbon=args.ribbon,
            reg_name=args.reg_name,
            nobandingcorr=args.nobandingcorr,
            outdir=args.outdir,
            stages=stages,
            is_longitudinal=args.is_longitudinal,
            longitudinal_study_dir=longitudinal_study_dir,
            longitudinal_template=args.longitudinal_template,
            asl_params=asl_params,
        )
    except Exception as e:
        logging.error(f"Error processing subject {subject_dir}:\n {e}")
        raise e


if __name__ == "__main__":
    main()
