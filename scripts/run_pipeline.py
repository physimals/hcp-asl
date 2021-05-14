"""
This script performs the full minimal pre-processing ASL pipeline 
for the Human Connectome Project (HCP) ASL data.

This currently requires that the script is called followed by 
the directories of the subjects of interest and finally the 
name of the MT correction scaling factors image.
"""

import sys
import os
from itertools import product
import logging

from hcpasl import __version__, __timestamp__, __sha1__
from hcpasl.distortion_correction import gradunwarp_and_topup
from hcpasl.m0_correction import correct_M0
from hcpasl.asl_correction import single_step_resample_to_asl0, single_step_resample_to_aslt1w
from hcpasl.asl_differencing import tag_control_differencing
from hcpasl.utils import setup_logger, create_dirs, split_mbpcasl
from hcpasl.qc import create_qc_report
from hcpasl.key_outputs import copy_key_outputs
from pathlib import Path
import subprocess
import argparse
from multiprocessing import cpu_count
import nibabel as nb

def process_subject(studydir, subid, mt_factors, mbpcasl, structural, 
                    fmaps, gradients, wmparc, ribbon, wbdir, use_t1=False, 
                    pvcorr=False, cores=cpu_count(), interpolation=3,
                    nobandingcorr=False, outdir="hcp_asl"):
    """
    Run the hcp-asl pipeline for a given subject.

    Parameters
    ----------
    studydir : pathlib.Path
        Path to the study's base directory.
    subid : str
        Subject id for the subject of interest.
    mt_factors : pathlib.Path
        Path to a .txt file of pre-calculated MT correction 
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
    wbdir : str
        path to development version of wb_command's bin directory 
        e.g. workbench/bin_macosx64
    use_t1 : bool, optional
        Whether or not to use the estimated T1 map in the 
        oxford_asl run in structural space.
    pvcorr : bool, optional
        Whether or not to run oxford_asl using pvcorr when 
        performing perfusion estimation in (ASL-gridded) T1 
        space.
    cores : int, optional
        Number of cores to use.
        When applying motion correction, this is the number 
        of cores that will be used by regtricks. Default is 
        the number of cores on your machine.
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
        Name of the main results directory. Default is 'hcp_asl'.
    """

    subject_dir = (studydir / subid).resolve(strict=True)
    logger = logging.getLogger("HCPASL")

    # initial set-up for the pipeline: create results directories
    logger.info("Creating main results directories.")
    asl_dir, aslt1w_dir = [subject_dir/outdir/name for name in ("ASL", "T1w/ASL")]
    tis_dir, calib0_dir, calib1_dir = [asl_dir/name for name in ("TIs", "Calib/Calib0", "Calib/Calib1")]
    create_dirs([asl_dir, aslt1w_dir, tis_dir, calib0_dir, calib1_dir])
    # split mbPCASL sequence into TIs and calibration images
    logger.info("Splitting mbPCASL sequence into ASL series and calibration images.")
    tis_name, calib0_name, calib1_name = [d/name for d, name in zip((tis_dir, calib0_dir, calib1_dir), 
                                                                    ("tis.nii.gz", "calib0.nii.gz", "calib1.nii.gz"))]
    split_mbpcasl(mbpcasl, tis_name, calib0_name, calib1_name)

    # run gradient_unwarp and topup, storing results 
    # in gradunwarp_dir and topup_dir respectively
    logger.info("Running gradient_unwarp and topup.")
    gradunwarp_dir = asl_dir/"gradient_unwarp"
    topup_dir = asl_dir/"topup"
    gradunwarp_and_topup(vol=str(calib0_name), 
                         coeffs_path=gradients, 
                         gradunwarp_dir=gradunwarp_dir, 
                         topup_dir=topup_dir, 
                         pa_sefm=str(fmaps["PA"]), 
                         ap_sefm=str(fmaps["AP"]), 
                         interpolation=interpolation)

    # apply corrections to the calibration images
    logger.info("Running M0 corrections.")
    hcppipedir = Path(os.environ["HCPPIPEDIR"])
    corticallut = hcppipedir/'global/config/FreeSurferCorticalLabelTableLut.txt'
    subcorticallut = hcppipedir/'global/config/FreeSurferSubcorticalLabelTableLut.txt'
    t1w_dir = structural["struct"].parent
    correct_M0(subject_dir=subject_dir, 
               calib_dir=calib0_dir.parent, 
               mt_factors=mt_factors, 
               t1w_dir=t1w_dir, 
               aslt1w_dir=aslt1w_dir, 
               gradunwarp_dir=gradunwarp_dir, 
               topup_dir=topup_dir, 
               wmparc=wmparc, 
               ribbon=ribbon, 
               corticallut=corticallut, 
               subcorticallut=subcorticallut, 
               interpolation=interpolation, 
               nobandingcorr=nobandingcorr, 
               outdir=outdir)

    # correct ASL series for distortion, bias, motion and banding
    # giving an ASL series in ASL0 space
    logger.info("Estimating ASL motion.")
    bias_field = calib0_dir/"BiasCorr/calib0_bias.nii.gz"
    if not nobandingcorr:
        calib_corr = calib0_dir/"MTCorr/calib0_mtcorr_gdc_dc_restore.nii.gz"
    else:
        calib_corr = calib0_dir/"BiasCorr/calib0_gdc_dc_restore.nii.gz"
    calib2struct = calib0_dir/"DistCorr/asl2struct.mat"
    single_step_resample_to_asl0(subject_dir=subject_dir, 
                                 tis_dir=tis_dir, 
                                 mt_factors=mt_factors, 
                                 bias_name=bias_field, 
                                 calib_name=calib_corr, 
                                 calib2struct=calib2struct, 
                                 gradunwarp_dir=gradunwarp_dir, 
                                 topup_dir=topup_dir, 
                                 t1w_dir=t1w_dir, 
                                 cores=cores, 
                                 interpolation=interpolation, 
                                 nobandingcorr=nobandingcorr, 
                                 outdir=outdir)
                                 
    # perform tag-control subtraction in ASL0 space
    logger.info("Performing tag-control subtraction of the corrected ASL series in ASL0 space.")
    if not nobandingcorr:
        series = tis_dir/"tis_gdc_dc_moco_restore_bandcorr.nii.gz"
    else:
        series = tis_dir/"tis_gdc_dc_moco_restore.nii.gz"
    scaling_factors = tis_dir/"combined_scaling_factors.nii.gz"
    betas_dir = tis_dir/"Betas"
    tag_control_differencing(series, scaling_factors, betas_dir, subject_dir, outdir)

    # estimate perfusion in ASL0 space using oxford_asl
    logger.info("Running oxford_asl in ASL0 space.")

    beta_perf = betas_dir/"beta_perf.nii.gz"
    asl0_brainmask = tis_dir/"aslfs_mask.nii.gz"
    oxford_asl_dir = tis_dir/"OxfordASL"
    oxford_asl_dir.mkdir(exist_ok=True)
    logger_oxasl = setup_logger("HCPASL.oxford_asl", oxford_asl_dir/"oxford_asl.log", "INFO")
    oxford_asl_call = [
        "oxford_asl",
        f"-i {str(betas_dir/'beta_perf.nii.gz')}", f"-o {str(oxford_asl_dir)}",
        f"-m {str(asl0_brainmask)}", "--tis=1.7,2.2,2.7,3.2,3.7", 
        "--slicedt=0.059", "--sliceband=10", "--casl", 
        "--ibf=tis", "--iaf=diff", "--rpts=6,6,6,10,15",
        "--fixbolus", "--bolus=1.5", "--te=19",
        "--debug", "--spatial=off"
    ]
    if use_t1:
        est_t1 = tis_dir/"SatRecov2/spatial/mean_T1t_filt.nii.gz"
        oxford_asl_call.append(f"--t1im={str(est_t1)}")
    oxford_asl_call = " ".join(oxford_asl_call)
    logger_oxasl.info(oxford_asl_call)
    process = subprocess.Popen(oxford_asl_call, shell=True, stdout=subprocess.PIPE)
    while 1:
        retcode = process.poll()
        line = process.stdout.readline().decode("utf-8")
        logger.info(line)
        if line == "" and retcode is not None:
            break
    if retcode != 0:
        logger.info(f"retcode={retcode}")
        logger.exception("Process failed.")
        
    # get data in ASLT1w space
    logger.info("Get data into ASLT1w space and re-estimate bias field.")
    if not nobandingcorr:
        asl_scaling_factors = tis_dir/"STCorr2/combined_scaling_factors_asln.nii.gz"
        mt_name = mt_factors
    else:
        asl_scaling_factors, mt_name = None, None
    t1_est = tis_dir/"SatRecov2/spatial/mean_T1t_filt.nii.gz"
    single_step_resample_to_aslt1w(asl_name=tis_name,
                                   calib_name=calib0_name,
                                   subject_dir=subject_dir,
                                   t1w_dir=t1w_dir,
                                   aslt1w_dir=aslt1w_dir,
                                   moco_dir=tis_dir/"MoCo/asln2m0_final.mat",
                                   perfusion_name=tis_dir/"OxfordASL/native_space/perfusion.nii.gz",
                                   gradunwarp_dir=gradunwarp_dir,
                                   topup_dir=topup_dir,
                                   ribbon=ribbon,
                                   wmparc=wmparc,
                                   corticallut=corticallut,
                                   subcorticallut=subcorticallut,
                                   asl_scaling_factors=asl_scaling_factors,
                                   mt_factors=mt_name,
                                   t1_est=t1_est,
                                   nobandingcorr=nobandingcorr,
                                   interpolation=interpolation,
                                   cores=cores)

    # perform partial volume estimation
    logger.info("Performing partial volume estimation.")
    pves_dir = aslt1w_dir/"PVEs"
    pves_dir.mkdir(exist_ok=True)
    logger_pv = setup_logger("HCPASL.pv_est", pves_dir/"pv_est.log", "INFO")
    pv_est_call = [
        "pv_est",
        str(subject_dir.parent),
        subject_dir.stem,
        "--cores", str(cores),
        "--outdir", outdir,
        "--interpolation", str(interpolation)
    ]
    process = subprocess.Popen(pv_est_call, stdout=subprocess.PIPE)
    while 1:
        retcode = process.poll()
        line = process.stdout.readline().decode("utf-8")
        logger.info(line)
        if line == "" and retcode is not None:
            break
    if retcode != 0:
        logger.info(f"retcode={retcode}")
        logger.exception("Process failed.")

    
    # perform tag-control subtraction in ASLT1w space
    logger.info("Performing tag-control subtraction of the corrected ASL series in ASLT1w space.")
    aslt1w_dir = aslt1w_dir
    series = aslt1w_dir/"TIs/asl_corr.nii.gz"
    scaling_factors = aslt1w_dir/"TIs/combined_scaling_factors.nii.gz"
    betas_dir = aslt1w_dir/"TIs/Betas"
    tag_control_differencing(series, scaling_factors, betas_dir, subject_dir, outdir)

    # final perfusion estimation in ASLT1w space
    logger.info("Running oxford_asl in ASLT1w space.")
    pve_dir = aslt1w_dir/"PVEs"
    oxford_aslt1w_dir = aslt1w_dir/"TIs/OxfordASL"
    oxford_aslt1w_dir.mkdir(exist_ok=True)
    logger_oxaslt1w = setup_logger("HCPASL.oxford_aslt1w", oxford_aslt1w_dir/"oxford_aslt1w.log", "INFO")
    oxford_aslt1w_call = [
        "oxford_asl",
        f"-i {str(betas_dir/'beta_perf.nii.gz')}",
        f"-o {str(oxford_aslt1w_dir)}",
        f"--pvgm={str(pve_dir/'pve_GM.nii.gz')}",
        f"--pvwm={str(pve_dir/'pve_WM.nii.gz')}",
        f"--csf={str(pve_dir/'vent_csf_mask.nii.gz')}",
        f"-c {str(aslt1w_dir/'Calib/Calib0/calib0_corr_aslt1w.nii.gz')}",
        f"-m {str(aslt1w_dir/'TIs/reg/ASL_FoV_brain_mask.nii.gz')}",
        f"--tiimg={str(aslt1w_dir/'TIs/timing_img_aslt1w.nii.gz')}",
        "--casl",
        "--ibf=tis",
        "--iaf=diff",
        "--rpts=6,6,6,10,15",
        "--fixbolus",
        "--bolus=1.5",
        "--te=19",
        "--debug",
        "--spatial=off",
        "--tr=8"
    ]
    if pvcorr:
        oxford_aslt1w_call.append("--pvcorr")
    if use_t1:
        est_t1 = aslt1w_dir/"TIs/reg/mean_T1t_filt_aslt1w.nii.gz"
        oxford_aslt1w_call.append(f"--t1im={str(est_t1)}")
    oxford_aslt1w_call = " ".join(oxford_aslt1w_call)
    logger_oxaslt1w.info(oxford_aslt1w_call)
    process = subprocess.Popen(oxford_aslt1w_call, shell=True, stdout=subprocess.PIPE)
    while 1:
        retcode = process.poll()
        line = process.stdout.readline().decode("utf-8")
        logger.info(line)
        if line == "" and retcode is not None:
            break
    if retcode != 0:
        logger.info(f"retcode={retcode}")
        logger.exception("Process failed.")

    logger.info("Projecting volumetric results to surface.")
    project_to_surface(studydir, subid, outdir=outdir, wbdir=wbdir)

    logger.info("Copying key outputs to $StudyDir/$SubID/T1w/ASL and $StudyDir/$SubID/MNINonLinear/ASL")
    copy_outputs(studydir, subid, outdir)

    logger.info("Creating QC report.")
    create_qc_report(subject_dir, outdir)

def project_to_surface(studydir, subid, outdir, wbdir, lowresmesh="32", FinalASLRes="2.5", 
                       SmoothingFWHM="2", GreyOrdsRes="2", RegName="MSMSulc"):
    """
    Project perfusion results to the cortical surface and generate
    CIFTI representation which includes both low res mesh surfaces
    in MSMSulc Atlas space, and subcortical structures in MNI 
    voxel space

    Parameters
    ----------
    studydir : pathlib.Path
        Path to the study's base directory.
    subid : str
        Subject id for the subject of interest.
    """
    logger = setup_logger("HCPASL.project", studydir/subid/outdir/"project.log", "INFO")
    # Projection scripts path:
    script         = "PerfusionCIFTIProcessingPipeline.sh"
    wb_path        = str(Path(wbdir).resolve(strict=True))

    ASLVariable    = ["perfusion_calib", "arrival"]
    ASLVariableVar = ["perfusion_var_calib", "arrival_var"]

    for idx in range(2):
        non_pvcorr_cmd = [script, studydir, subid, ASLVariable[idx], ASLVariableVar[idx], lowresmesh,
                FinalASLRes, SmoothingFWHM, GreyOrdsRes, RegName, wb_path, "false", outdir]

        pvcorr_cmd = [script, studydir, subid, ASLVariable[idx], ASLVariableVar[idx], lowresmesh,
                FinalASLRes, SmoothingFWHM, GreyOrdsRes, RegName, wb_path, "true", outdir]
        
        process = subprocess.Popen(non_pvcorr_cmd, stdout=subprocess.PIPE)
        while 1:
            retcode = process.poll()
            line = process.stdout.readline().decode("utf-8")
            logger.info(line)
            if line == "" and retcode is not None:
                break
        if retcode != 0:
            logger.info(f"retcode={retcode}")
            logger.exception("Process failed.")
        process = subprocess.Popen(pvcorr_cmd, stdout=subprocess.PIPE)
        while 1:
            retcode = process.poll()
            line = process.stdout.readline().decode("utf-8")
            logger.info(line)
            if line == "" and retcode is not None:
                break
        if retcode != 0:
            logger.info(f"retcode={retcode}")
            logger.exception("Process failed.")

def copy_outputs(studydir, subid, outdir):
    """
    Copy key pipeline outputs to the T1w and MNI aligned high level ASL directory

    Parameters
    ----------
    studydir : pathlib.Path
        Path to the study's base directory.
    subid : str
        Subject id for the subject of interest.
    """
       
    path_to_outs   = str(studydir/subid/outdir)
    mni_raw = str(studydir/subid/f"{subid}_V1_MR/resources/Structural_preproc/files/{subid}_V1_MR/MNINonLinear")
    t1w_preproc = str(studydir/subid/f"{subid}_V1_MR/resources/Structural_preproc/files/{subid}_V1_MR/T1w")
    copy_key_outputs(path_to_outs, t1w_preproc, mni_raw)




def main():
    """
    Main entry point for the hcp-asl pipeline.
    """
    # argument handling
    parser = argparse.ArgumentParser(
        description="This script performs the minimal processing for the "
                    + "HCP-Aging ASL data.")
    parser.add_argument(
        "--studydir",
        help="Path to the study's base directory.",
        required=True
    )
    parser.add_argument(
        "--subid",
        help="Subject id for the subject of interest.",
        required=True
    )
    parser.add_argument(
        "--mtname",
        help="Filename of the empirically estimated MT-correction"
            + "scaling factors.",
        required=not "--nobandingcorr" in sys.argv
    )
    parser.add_argument(
        "-g",
        "--grads",
        help="Filename of the gradient coefficients for gradient"
            + "distortion correction.",
        required=True
    )
    parser.add_argument(
        "-s",
        "--struct",
        help="Filename for the acpc-aligned, dc-restored structural image.",
        required=True
    )
    parser.add_argument(
        "--sbrain",
        help="Filename for the brain-extracted acpc-aligned, "
            + "dc-restored structural image.",
        required=True
    )
    parser.add_argument(
        "--mbpcasl",
        help="Filename for the mbPCASLhr acquisition.",
        required=True
    )
    parser.add_argument(
        "--fmap_ap",
        help="Filename for the AP fieldmap for use in distortion correction",
        required=True
    )
    parser.add_argument(
        "--fmap_pa",
        help="Filename for the PA fieldmap for use in distortion correction",
        required=True
    )
    parser.add_argument(
        '--use_t1',
        help="If this flag is provided, the T1 estimates from the satrecov "
            + "will also be registered to ASL-gridded T1 space for use in "
            + "perfusion estimation via oxford_asl.",
        action='store_true'
    )
    parser.add_argument(
        '--pvcorr',
        help="If this flag is provided, oxford_asl will be run using the "
            + "--pvcorr flag.",
        action='store_true'
    )
    parser.add_argument(
        '--wmparc',
        help="wmparc.nii.gz from FreeSurfer for use in SE-based bias correction.",
        default=None,
        required=True
    )
    parser.add_argument(
        '--ribbon',
        help="ribbon.nii.gz from FreeSurfer for use in SE-based bias correction.",
        default=None,
        required=True
    )
    parser.add_argument(
        "-c",
        "--cores",
        help="Number of cores to use when applying motion correction and "
            +"other potentially multi-core operations. Default is the "
            +f"number of cores your machine has ({cpu_count()}).",
        default=cpu_count(),
        type=int,
        choices=range(1, cpu_count()+1)
    )
    parser.add_argument(
        "--interpolation",
        help="Interpolation order for registrations. This can be any "
            +"integer from 0-5 inclusive. Default is 3. See scipy's "
            +"map_coordinates for more details.",
        default=3,
        type=int,
        choices=range(0, 5+1)
    )
    parser.add_argument(
        "--nobandingcorr",
        help="If this option is provided, the MT and ST banding corrections "
            +"won't be applied. This is to be used to compare the difference "
            +"our banding corrections make.",
        action="store_true"
    )
    parser.add_argument(
        "--fabberdir",
        help="User Fabber executable in <fabberdir>/bin/ for users"
            + "with FSL < 6.0.4"
    )
    parser.add_argument(
        "--wbdir",
        help="Location of wb_command/bin_macosx64 (>= v1.5.0).",
        required=True
    )
    parser.add_argument(
        "--outdir",
        help="Name of the directory within which we will store all of the "
            +"pipeline's outputs in sub-directories. Default is the subject's"
            +"base directory.",
        default=""
    )
    parser.add_argument(
        "-v", "--verbose",
        help="If this option is provided, stdout will go to the terminal "
            +"as well as to a logfile. Default is False.",
        action="store_true"
    )
    # assign arguments to variables
    args = parser.parse_args()
    studydir = Path(args.studydir).resolve(strict=True)
    subid = args.subid

    # set up logging
    # create file handler
    base_dir = Path(studydir/subid/args.outdir)
    base_dir.mkdir(exist_ok=True)
    fh_name = base_dir/f"{subid}.log"
    logger = setup_logger("HCPASL", fh_name, "INFO", args.verbose)
    
    logger.info(f"Welcome to HCPASL v{__version__} (commit {__sha1__} on {__timestamp__}.")
    logger.info(args)

    # parse remaining arguments
    if args.mtname:
        mtname = Path(args.mtname).resolve(strict=True)
    else:
        mtname = None
    structural = {'struct': Path(args.struct).resolve(strict=True),
                  'sbrain': Path(args.sbrain).resolve(strict=True)}
    mbpcasl = Path(args.mbpcasl).resolve(strict=True)
    fmaps = {
        'AP': Path(args.fmap_ap).resolve(strict=True), 
        'PA': Path(args.fmap_pa).resolve(strict=True)
    }
    grads = Path(args.grads).resolve(strict=True)

    if args.fabberdir:
        if not os.path.isfile(os.path.join(args.fabberdir, "bin", "fabber_asl")):
            logger.error("ERROR: specified Fabber in %s, but no fabber_asl executable found in %s/bin" % (args.fabberdir, args.fabberdir))
            sys.exit(1)

        # To use a custom Fabber executable we set the FSLDEVDIR environment variable
        # which prioritises executables in $FSLDEVDIR/bin over those in $FSLDIR/bin.
        # Note that this could cause problems in the unlikely event that the user
        # already has a $FSLDEVDIR set up with custom copies of other things that
        # oxford_asl uses...
        logger.info("Using Fabber-ASL executable %s/bin/fabber_asl" % args.fabberdir)
        os.environ["FSLDEVDIR"] = os.path.abspath(args.fabberdir)

    # process subject
    logger.info(f"Processing subject {studydir/subid}.")
    process_subject(studydir=studydir,
                    subid=subid,
                    mt_factors=mtname,
                    cores=args.cores,
                    interpolation=args.interpolation,
                    gradients=grads,
                    mbpcasl=mbpcasl,
                    structural=structural,
                    fmaps=fmaps,
                    use_t1=args.use_t1,
                    pvcorr=args.pvcorr,
                    wmparc=args.wmparc,
                    ribbon=args.ribbon,
                    nobandingcorr=args.nobandingcorr,
                    outdir=args.outdir,
                    wbdir=args.wbdir
                    )

if __name__ == '__main__':
    main()