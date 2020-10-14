"""
This script performs the full minimal pre-processing ASL pipeline 
for the Human Connectome Project (HCP) ASL data.

This currently requires that the script is called followed by 
the directories of the subjects of interest and finally the 
name of the MT correction scaling factors image.
"""

import sys
import os

from hcpasl.initial_bookkeeping import initial_processing
from hcpasl.m0_mt_correction import correct_M0
from hcpasl.asl_correction import hcp_asl_moco
from hcpasl.asl_differencing import tag_control_differencing
from hcpasl.asl_perfusion import run_fabber_asl, run_oxford_asl
from hcpasl.projection import project_to_surface
from pathlib import Path
import subprocess
import argparse
from multiprocessing import cpu_count

def process_subject(subject_dir, mt_factors, cores, order, mbpcasl, structural, surfaces, 
                    fmaps, gradients=None, use_t1=False, pvcorr=False, use_sebased=False,
                    wmparc=None, ribbon=None, cLUT=None, scLUT=None):
    """
    Run the hcp-asl pipeline for a given subject.

    Parameters
    ----------
    subject_dir : str
        Path to the subject's base directory.
    mt_factors : str
        Path to a .txt file of pre-calculated MT correction 
        factors.
    cores : int
        Number of cores to use.
        When applying motion correction, this is the number 
        of cores that will be used by regtricks.
    order : int
        The interpolation order to use for registrations.
        Regtricks passes this on to scipy's map_coordinates. 
        The meaning of the value can be found in the scipy 
        documentation.
    mbpcasl : str
        Path to the subject's mbPCASL sequence.
    structural : dict
        Contains the locations of important structural files.
    surfaces : dict
        Contains the locations of the surfaces needed for the 
        pipeline.
    fmaps : dict
        Contains the locations of the fieldmaps needed for 
        distortion correction.
    gradients : str, optional
        Path to a gradient coefficients file for use in 
        gradient distortion correction.
    use_t1 : bool, optional
        Whether or not to use the estimated T1 map in the 
        oxford_asl run in structural space.
    use_t1 : bool, optional
        Whether or not to run oxford_asl using pvcorr when 
        performing perfusion estimation in (ASL-gridded) T1 
        space.
    use_sebased : bool, optional
        Whether or not to use HCP's SE-based bias-correction 
        to refine the bias correction obtained using FAST at 
        the beginning of the pipeline. Default is False.
    """
    subject_dir = Path(subject_dir)
    mt_factors = Path(mt_factors)
    initial_processing(subject_dir, mbpcasl=mbpcasl, structural=structural, surfaces=surfaces)
    correct_M0(subject_dir, mt_factors)
    hcp_asl_moco(subject_dir, mt_factors, cores=cores, order=order)
    for target in ('asl', 'structural'):
        dist_corr_call = [
            "hcp_asl_distcorr",
            str(subject_dir.parent),
            subject_dir.stem,
            "--target", target,
            "--fmap_ap", fmaps['AP'], 
            "--fmap_pa", fmaps['PA'],
            '--mt', mt_factors
        ]
        if gradients:
            dist_corr_call.append('--grads')
            dist_corr_call.append(gradients)
        if use_t1 and (target=='structural'):
            dist_corr_call.append('--use_t1')
        if use_sebased and (target=='structural'):
            dist_corr_call.append('--sebased')
        subprocess.run(dist_corr_call, check=True)
        if target == 'structural':
            pv_est_call = [
                "pv_est",
                str(subject_dir.parent),
                subject_dir.stem
            ]
            subprocess.run(pv_est_call, check=True)
        if use_sebased and (target=='structural'):
            calib_name = subject_dir/'T1w/ASL/Calib/Calib0/DistCorr/calib0_dcorr.nii.gz'
            asl_name = subject_dir/'T1w/ASL/TIs/DistCorr/tis_distcorr.nii.gz'
            mask_name = subject_dir/'T1w/ASL/reg/ASL_grid_T1w_acpc_dc_restore_brain_mask.nii.gz'
            fmapmag_name = subject_dir/'T1w/ASL/reg/fmap/fmapmag_aslstruct.nii.gz'
            out_dir = subject_dir/'T1w/ASL/TIs/BiasCorr'
            sebased_cmd = [
                'get_sebased_bias',
                '-i', calib_name,
                '--asl', asl_name,
                '-f', fmapmag_name,
                '-m', mask_name,
                '--wmparc', wmparc,
                '--ribbon', ribbon,
                '--cortical', cLUT,
                '--subcortical', scLUT,
                '-o', out_dir,
                '--debug'
            ]
            subprocess.run(sebased_cmd, check=True)
        if use_sebased and (target=='structural'):
            series = subject_dir/'T1w/ASL/TIs/BiasCorr/tis_secorr.nii.gz'
        elif target=='structural':
            series = subject_dir/'T1w/ASL/TIs/DistCorr/tis_distcorr.nii.gz'
        else:
            series = subject_dir/'ASL/TIs/DistCorr/tis_distcorr.nii.gz'
        tag_control_differencing(series, subject_dir, target=target)
        run_oxford_asl(subject_dir, target=target, use_t1=use_t1, pvcorr=pvcorr)
        project_to_surface(subject_dir, target=target)

def main():
    """
    Main entry point for the hcp-asl pipeline.
    """
    # argument handling
    parser = argparse.ArgumentParser(
        description="This script performs the minimal processing for the "
                    + "HCP-Aging ASL data.")
    parser.add_argument(
        "subject_dir",
        help="The directory of the subject you wish to process."
    )
    parser.add_argument(
        "scaling_factors",
        help="Filename of the empirically estimated MT-correction"
            + "scaling factors."
    )
    parser.add_argument(
        "-g",
        "--grads",
        help="Filename of the gradient coefficients for gradient"
            + "distortion correction (optional)."
    )
    parser.add_argument(
        "-s",
        "--struct",
        help="Filename for the acpc-aligned, dc-restored structural image."
    )
    parser.add_argument(
        "--sbrain",
        help="Filename for the brain-extracted acpc-aligned, "
            + "dc-restored structural image."
    )
    parser.add_argument(
        "--lmid",
        help="Filename for the left mid surface."
    )
    parser.add_argument(
        "--rmid",
        help="Filename for the right mid surface."
    )
    parser.add_argument(
        "--lwhite",
        help="Filename for the left white surface."
    )
    parser.add_argument(
        "--rwhite",
        help="Filename for the right white surface."
    )
    parser.add_argument(
        "--lpial",
        help="Filename for the left pial surface."
    )
    parser.add_argument(
        "--rpial",
        help="Filename for the right pial surface."
    )
    parser.add_argument(
        "-i",
        "--input",
        help="Filename for the mbPCASLhr acquisition."
    )
    parser.add_argument(
        "--fmap_ap",
        help="Filename for the AP fieldmap for use in distortion correction"
    )
    parser.add_argument(
        "--fmap_pa",
        help="Filename for the PA fieldmap for use in distortion correction"
    )
    parser.add_argument(
        '--use_t1',
        help="If this flag is provided, the T1 estimates from the satrecov "
            + "will also be registered to ASL-gridded T1 space for use in "
            + "perfusion estimation via oxford_asl.",
        action='store_true'
    )
    parser.add_argument(
        '--sebased',
        help="If this flag is provided, the distortion warps and motion "
            +"estimates will be applied to the MT-corrected but not bias-"
            +"corrected calibration and ASL images. The bias-field will "
            +"then be estimated from the calibration image using HCP's "
            +"SE-based algorithm and applied in subsequent steps.",
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
        help="wmparc.mgz from FreeSurfer",
        default=None
    )
    parser.add_argument(
        '--ribbon',
        help="ribbon.mgz from FreeSurfer",
        default=None
    )
    parser.add_argument(
        "--cortical",
        help="Filename for FreeSurfer's Cortical Lable Table",
        default=None
    )
    parser.add_argument(
        "--subcortical",
        help="Filename for FreeSurfer's Subcortical Lable Table",
        default=None
    )
    parser.add_argument(
        "-c",
        "--cores",
        help="Number of cores to use for registration operations. "
            + f"Your PC has {cpu_count()}. Default is 1.",
        default=1,
        type=int
    )
    parser.add_argument(
        "--interpolation",
        help="Interpolation order for registrations. Default is 3.",
        default=3,
        type=int
    )
    parser.add_argument(
        "--fabberdir",
        help="User Fabber executable in <fabberdir>/bin/ for users"
            + "with FSL < 6.0.4"
    )
    # assign arguments to variables
    args = parser.parse_args()
    mt_name = args.scaling_factors
    subject_dir = args.subject_dir
    structural = {'struct': args.struct, 'sbrain': args.sbrain}
    surfaces = {
        'L_mid': args.lmid, 'R_mid': args.rmid,
        'L_white': args.lwhite, 'R_white':args.rwhite,
        'L_pial': args.lpial, 'R_pial': args.rpial
    }
    mbpcasl = args.input
    fmaps = {'AP': args.fmap_ap, 'PA': args.fmap_pa}
    use_t1 = args.use_t1
    use_sebased = args.sebased
    pvcorr = args.pvcorr
    wmparc = args.wmparc
    ribbon = args.ribbon
    scLUT = args.subcortical
    cLUT = args.cortical
    cores = args.cores
    order = args.interpolation
    if args.fabberdir:
        if not os.path.isfile(os.path.join(args.fabberdir, "bin", "fabber_asl")):
            print("ERROR: specified Fabber in %s, but no fabber_asl executable found in %s/bin" % (args.fabberdir, args.fabberdir))
            sys.exit(1)

        # To use a custom Fabber executable we set the FSLDEVDIR environment variable
        # which prioritises executables in $FSLDEVDIR/bin over those in $FSLDIR/bin.
        # Note that this could cause problems in the unlikely event that the user
        # already has a $FSLDEVDIR set up with custom copies of other things that
        # oxford_asl uses...
        print("Using Fabber-ASL executable %s/bin/fabber_asl" % args.fabberdir)
        os.environ["FSLDEVDIR"] = os.path.abspath(args.fabberdir)

    print(f"Processing subject {subject_dir}.")
    if args.grads:
        print("Including gradient distortion correction step.")
        process_subject(subject_dir=subject_dir,
                        mt_factors=mt_name,
                        cores=cores,
                        order=order,
                        gradients=args.grads,
                        mbpcasl=mbpcasl,
                        structural=structural,
                        surfaces=surfaces,
                        fmaps=fmaps,
                        use_t1=use_t1,
                        pvcorr=pvcorr,
                        use_sebased=use_sebased,
                        wmparc=wmparc,
                        ribbon=ribbon,
                        cLUT=cLUT,
                        scLUT=scLUT
                        )
    else:
        print("Not including gradient distortion correction step.")
        process_subject(subject_dir=subject_dir,
                        mt_factors=mt_name,
                        cores=cores,
                        order=order,
                        mbpcasl=mbpcasl,
                        structural=structural,
                        surfaces=surfaces,
                        fmaps=fmaps,
                        use_t1=use_t1,
                        pvcorr=pvcorr,
                        use_sebased=use_sebased,
                        wmparc=wmparc,
                        ribbon=ribbon,
                        cLUT=cLUT,
                        scLUT=scLUT
                        )

if __name__ == '__main__':
    main()