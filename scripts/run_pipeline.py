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
from shutil import copy, rmtree

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
    territories_atlas,
    territories_labels,
    use_t1=False,
    cores=1,
    interpolation=3,
    nobandingcorr=False,
    outdir="hcp_asl",
    stages=set(range(14)),
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
    """

    if not isinstance(stages, (list, set)):
        raise RuntimeError("stages must either be a set or list of ints")

    # create results directories
    logging.info("Creating main results directories.")
    asl_dir, aslt1w_dir = [subject_dir / outdir / name for name in ("ASL", "T1w/ASL")]
    label_control_dir, calib0_dir, calib1_dir = [
        asl_dir / name
        for name in ("label_control", "calibration/calib0", "calibration/calib1")
    ]
    for d in [asl_dir, aslt1w_dir, label_control_dir, calib0_dir, calib1_dir]:
        d.mkdir(exist_ok=True, parents=True)

    # split ASL sequence into label-control label_control and calibration images
    tis_name, calib0_name, calib1_name = [
        d / name
        for d, name in zip(
            (label_control_dir, calib0_dir, calib1_dir),
            ("label_control.nii.gz", "calib0.nii.gz", "calib1.nii.gz"),
        )
    ]
    if 0 in stages:
        logging.info(
            "Stage 0: Splitting ASL sequence into label-control pairs and calibration images."
        )
        split_asl(mbpcasl, tis_name, calib0_name, calib1_name)

    # Run gradient_unwarp and topup on the calibration images
    # results stored in gradunwarp_dir and topup_dir respectively
    gradunwarp_dir = asl_dir / "gradient_unwarp"
    topup_dir = asl_dir / "topup"
    gd_corr = gradients is not None
    if 1 in stages:
        logging.info(
            "Stage 1: Derive gradient and susceptibility distortion correction."
        )
        if not gd_corr:
            logging.info(
                "Gradient coefficient file not provided, derivation will be skipped."
            )
            logging.info("Deriving susceptibility distortion correction only.")
        derive_gdc_sdc(
            vol=str(calib0_name),
            coeffs_path=gradients,
            gradunwarp_dir=gradunwarp_dir,
            topup_dir=topup_dir,
            pa_sefm=str(fmaps["PA"]),
            ap_sefm=str(fmaps["AP"]),
            interpolation=interpolation,
            gd_corr=gd_corr,
        )

    # Apply corrections derived thus far to M0 image
    # Estimate biasfield
    t1w_dir = structural["struct"].parent
    if 2 in stages:
        logging.info("Stage 2: Derive and apply initial corrections to M0 image.")
        initial_corrections_calibration(
            subject_id=subid,
            calib_dir=calib0_dir.parent,
            eb_factors=eb_factors,
            t1w_dir=t1w_dir,
            aslt1w_dir=aslt1w_dir,
            gradunwarp_dir=gradunwarp_dir,
            gd_corr=gd_corr,
            topup_dir=topup_dir,
            wmparc=wmparc,
            ribbon=ribbon,
            interpolation=interpolation,
            nobandingcorr=nobandingcorr,
        )

    calib_corr = calib0_dir / "calib0_initial_corrected.nii.gz"
    bias_field = calib0_dir / "bias_correction/calib0_biasfield.nii.gz"
    calib2struct = calib0_dir / "registration/asl2struct.mat"

    # Apply corrections derived thus far to ASL timeseries
    if 3 in stages:
        logging.info("Stage 3: Derive and apply initial corrections to ASL timeseries.")
        initial_corrections_asl(
            subject_dir=subject_dir,
            label_control_dir=label_control_dir,
            eb_factors=eb_factors,
            bias_name=bias_field,
            calib_name=calib_corr,
            calib2struct=calib2struct,
            gradunwarp_dir=gradunwarp_dir,
            gd_corr=gd_corr,
            topup_dir=topup_dir,
            t1w_dir=t1w_dir,
            cores=cores,
            interpolation=interpolation,
            nobandingcorr=nobandingcorr,
        )

    # perform tag-control subtraction in ASL0 space
    asl_lc = label_control_dir / "label_control_corrected.nii.gz"
    scaling_factors = label_control_dir / "combined_scaling_factors_mc.nii.gz"
    asl_subtract = label_control_dir / "motion_subtraction"
    asl0_brainmask = label_control_dir / "brain_fov_mask.nii.gz"
    if 4 in stages:
        logging.info("Stage 4: Label-control subtraction in native ASL space.")
        tag_control_differencing(
            asl_lc, scaling_factors, asl_subtract, mask=asl0_brainmask
        )

    # estimate perfusion in ASL0 space using oxford_asl
    oxford_asl_dir = label_control_dir.parent / "perfusion_estimation"
    oxford_asl_dir.mkdir(exist_ok=True, parents=True)
    if 5 in stages:
        logging.info("Stage 5: Perfusion estimation in ASL native space.")
        logging.info(
            f"Copying oxford_asl inputs to one location ({str(oxford_asl_dir / 'oxford_asl_inputs')})."
        )
        oxasl_inputs = {
            "-i": asl_subtract / "beta_perf.nii.gz",
            "-m": asl0_brainmask,
        }
        copy_oxford_asl_inputs(oxasl_inputs, oxford_asl_dir / "oxford_asl_inputs")
        oxford_asl_call = [
            "oxford_asl",
            *[f"{k}={str(v)}" for k, v in oxasl_inputs.items()],
            f"-o={str(oxford_asl_dir)}",
            "--tis=" + ",".join([str(t) for t in TIS]),
            f"--slicedt={SLICEDT}",
            f"--sliceband={SLICEBAND}",
            "--casl",
            f"--ibf={IBF}",
            "--iaf=diff",
            "--rpts=" + ",".join([str(r) for r in RPTS]),
            "--fixbolus",
            f"--bolus={BOLUS}",
            f"--te={TE}",
            "--spatial=off",
            "--debug",
        ]
        if use_t1:
            est_t1 = (
                label_control_dir
                / "saturation_recovery/second/spatial/mean_T1t_filt.nii.gz"
            )
            oxford_asl_call.append(f"--t1im={str(est_t1)}")
        logging.info(oxford_asl_call)
        sp_run(oxford_asl_call)

    # get data in ASLT1w space
    if not nobandingcorr:
        asl_scaling_factors = (
            label_control_dir
            / "slicetime_correction/second/combined_scaling_factors_asln.nii.gz"
        )
        mt_name = eb_factors
    else:
        asl_scaling_factors, mt_name = None, None
    t1_est = (
        label_control_dir / "saturation_recovery/second/spatial/mean_T1t_filt.nii.gz"
    )
    if 6 in stages:
        logging.info(
            "Stage 6: Fully-correct ASL and calibration into ASL-gridded T1w space."
        )
        fully_correct_asl_calibration_aslt1w(
            asl_name=tis_name,
            calib_name=calib0_name,
            subid=subid,
            subject_dir=subject_dir,
            t1w_dir=t1w_dir,
            aslt1w_dir=aslt1w_dir,
            moco_dir=label_control_dir / "motion_correction/asln2calibration_final.mat",
            perfusion_name=label_control_dir.parent
            / "perfusion_estimation/native_space/perfusion.nii.gz",
            gradunwarp_dir=gradunwarp_dir,
            gd_corr=gd_corr,
            topup_dir=topup_dir,
            ribbon=ribbon,
            wmparc=wmparc,
            asl_scaling_factors=asl_scaling_factors,
            eb_factors=mt_name,
            t1_est=t1_est,
            nobandingcorr=nobandingcorr,
            interpolation=interpolation,
            cores=cores,
        )
        copy(
            aslt1w_dir / "calibration/calib0/calib0_corrected.nii.gz",
            aslt1w_dir / "calib_corrected.nii.gz",
        )
        copy(
            aslt1w_dir / "label_control/label_control_corrected.nii.gz",
            aslt1w_dir / "label_control_corrected.nii.gz",
        )

    # perform partial volume estimation
    if 7 in stages:
        logging.info("Stage 7: Partial volume estimation in ASLT1w space.")
        run_pv_estimation(subject_dir, cores, outdir, interpolation)

    # perform tag-control subtraction in ASLT1w space
    aslt1w_dir = aslt1w_dir
    series = aslt1w_dir / "label_control/label_control_corrected.nii.gz"
    scaling_factors = aslt1w_dir / "label_control/label_control_scaling_factors.nii.gz"
    subtracted_dir = aslt1w_dir / "label_control/motion_subtraction"
    brainmask = aslt1w_dir / "registration/brain_fov_mask.nii.gz"
    if 8 in stages:
        logging.info("Stage 8: Label-control subtraction in ASLT1w space")
        tag_control_differencing(
            series, scaling_factors, subtracted_dir, mask=brainmask
        )
        copy(
            subtracted_dir / "beta_perf.nii.gz",
            aslt1w_dir / "label_control_corrected_subtracted.nii.gz",
        )

    # final perfusion estimation in ASLT1w space
    pve_dir = aslt1w_dir / "pvs"
    gm_pve, wm_pve = [pve_dir / f"pv_{tiss}.nii.gz" for tiss in ("GM", "WM")]
    oxford_aslt1w_dir = aslt1w_dir / "perfusion_estimation"
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
            "-c": aslt1w_dir / "calibration/calib0/calib0_corrected.nii.gz",
            "-m": aslt1w_dir / "registration/brain_fov_mask.nii.gz",
            "--tiimg": aslt1w_dir / "label_control/timing_img.nii.gz",
        }
        if use_t1:
            oxasl_inputs["--t1im"] = (
                aslt1w_dir / "registration/mean_T1t_filt_aslt1w.nii.gz"
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
            est_t1 = aslt1w_dir / "registration/mean_T1t_filt.nii.gz"
            oxford_aslt1w_call.append(f"--t1im={str(est_t1)}")
        sp_run(oxford_aslt1w_call)

    mninonlinear_name = subject_dir / "MNINonLinear"
    if 10 in stages:
        logging.info("Stage 10: Summary statistics within ROIs.")
        roi_stats(
            struct_name=structural["struct"],
            oxford_asl_dir=oxford_aslt1w_dir,
            gm_pve=gm_pve,
            wm_pve=wm_pve,
            std2struct_name=mninonlinear_name / "xfms/standard2acpc_dc.nii.gz",
            roi_stats_dir=aslt1w_dir / "roi_stats",
            territories_atlas=territories_atlas,
            territories_labels=territories_labels,
        )

    if 11 in stages:
        logging.info("Stage 11: Volume to surface projection.")
        surface_projection_stage(
            subject_dir=subject_dir, subject_id=subid, outdir=outdir
        )

    if 12 in stages:
        logging.info(
            "Stage 12: Copy key results into $outdir/T1w/ASL and $outdir/MNINonLinear/ASL"
        )
        copy_outputs(subject_dir, outdir)

    if 13 in stages:
        logging.info("Stage 13: Create QC workbench scene.")
        create_qc_report(subject_id=subid, subject_dir=subject_dir, outdir=outdir)


def surface_projection_stage(
    subject_id,
    subject_dir,
    outdir,
    lowresmesh="32",
    FinalASLRes="2.5",
    SmoothingFWHM="2",
    GreyOrdsRes="2",
    RegName="MSMAll",
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
            RegName,
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
            RegName,
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
        default=set(range(14)),
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

    logging.info("All pipeline arguments:")
    for k, v in vars(args).items():
        logging.info(f"{k}: {v}")

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
            nobandingcorr=args.nobandingcorr,
            outdir=args.outdir,
            stages=args.stages,
        )
    except Exception as e:
        logging.error(f"Error processing subject {subject_dir}:\n {e}")
        raise e


if __name__ == "__main__":
    main()
