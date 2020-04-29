from m0_mt_correction import load_json, update_json
from initial_bookkeeping import create_dirs
from pathlib import Path
import subprocess
import numpy as np

def run_oxford_asl(subject_dir):
    # load subject's json
    json_dict = load_json(subject_dir)

    # directory for oxford_asl results
    oxford_dir = Path(json_dict['structasl']) / 'TIs/OxfordASL'
    pvgm_name = Path(json_dict['structasl']) / 'PVEs/pve_GM.nii.gz'
    pvwm_name = Path(json_dict['structasl']) / 'PVEs/pve_WM.nii.gz'
    calib_name = Path(json_dict['structasl']) / 'Calib/Calib0/DistCorr/calib0_dcorr.nii.gz'
    cmd = [
        "oxford_asl",
        f"-i {json_dict['beta_perf']}",
        f"-o {str(oxford_dir)}",
        "--casl",
        "--ibf=tis",
        "--iaf=diff",
        "--tis=1.7,2.2,2.7,3.2,3.7",
        "--rpts=6,6,6,10,15",
        "--fixbolus",
        "--bolus=1.5",
        "--pvcorr",
        f"-c {str(calib_name)}",
        "--cmethod=single",
        f"--pvgm={str(pvgm_name)}",
        f"--pvwm={str(pvwm_name)}",
        "--te=19",
        "--debug",
        "--spatial=off",
        "--slicedt=0.059",
        "--sliceband=10",
        f"-s {json_dict['T1w_acpc']}",
        f"--sbrain={json_dict['T1w_acpc_brain']}"
    ]
    subprocess.run(" ".join(cmd), shell=True)

    # add oxford_asl directory to the json
    important_names = {
        "oxford_asl": str(oxford_dir)
    }
    update_json(important_names, json_dict)