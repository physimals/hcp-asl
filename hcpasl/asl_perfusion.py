from .m0_mt_correction import load_json, update_json
from .initial_bookkeeping import create_dirs
from pathlib import Path
import subprocess
import numpy as np

def run_oxford_asl(subject_dir, target='structural'):
    # load subject's json
    json_dict = load_json(subject_dir)

    # base oxford_asl options common to both cases
    cmd = [
        "oxford_asl",
        f"-i {json_dict['beta_perf']}",
        "--casl",
        "--ibf=tis",
        "--iaf=diff",
        "--tis=1.7,2.2,2.7,3.2,3.7",
        "--rpts=6,6,6,10,15",
        "--fixbolus",
        "--bolus=1.5",
        "--te=19",
        "--debug",
        "--spatial=off",
        "--slicedt=0.059",
        "--sliceband=10"
    ]
    # directory for oxford_asl results
    if target == 'asl':
        oxford_dir = Path(json_dict['TIs_dir']) / 'SecondPass/OxfordASL'
        brain_mask = Path(json_dict['structasl']) / 'reg/asl_vol1_mask_init.nii.gz'
        extra_args = [
            f"-o {str(oxford_dir)}",
            f"-m {str(brain_mask)}"
        ]
    else:
        structasl_dir = Path(json_dict['structasl'])
        oxford_dir = structasl_dir / 'TIs/OxfordASL'
        pvgm_name = structasl_dir / 'PVEs/pve_GM.nii.gz'
        pvwm_name = structasl_dir / 'PVEs/pve_WM.nii.gz'
        csf_mask_name = structasl_dir / 'PVEs/vent_csf_mask.nii.gz'
        calib_name = structasl_dir / 'Calib/Calib0/DistCorr/calib0_dcorr.nii.gz'
        brain_mask = structasl_dir / 'reg/ASL_grid_T1w_acpc_dc_restore_brain_mask.nii.gz'
        extra_args = [
            f"-o {str(oxford_dir)}",
            f"--pvgm={str(pvgm_name)}",
            f"--pvwm={str(pvwm_name)}",
            f"--csf={str(csf_mask_name)}",
            f"-c {str(calib_name)}",
            f"-m {str(brain_mask)}"
            ]
    for extra_arg in extra_args:
        cmd.append(extra_arg)
    print(" ".join(cmd))
    subprocess.run(" ".join(cmd), shell=True)

    # add oxford_asl directory to the json
    important_names = {
        "oxford_asl": str(oxford_dir)
    }
    update_json(important_names, json_dict)