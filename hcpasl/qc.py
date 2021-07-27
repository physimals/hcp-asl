import subprocess
from hcpasl.utils import setup_logger
from nbclient import execute
import nbformat
from nbparameterise import extract_parameters, parameter_values, replace_definitions

from pathlib import Path
import os

import regtricks as rt

from importlib.resources import path as resource_path
from . import resources

def create_qc_report(subject_dir, outdir):
    # get location of this file
    p = resource_path(resources, 'report_template.ipynb')
    with p as filename:
        template_nb = nbformat.read(filename, as_version=4)

        # extract original parameters from the template notebook
        orig_parameters = extract_parameters(template_nb)

        # replace with parameters for this subject and execute notebook
        new_parameters = parameter_values(orig_parameters,
                                        subject_dir=str(subject_dir),
                                        outdir=str(outdir))
        new_nb = replace_definitions(template_nb, new_parameters, execute=False)
        _ = execute(new_nb)

        # save notebook in subject's main output directory
        new_nb_name = Path(subject_dir)/outdir/"report.ipynb"
        with open(new_nb_name, "w") as f:
            nbformat.write(new_nb, f)

def roi_stats(struct_name, oxford_asl_dir, gm_pve, wm_pve, std2struct_name,
              roi_stats_dir, territories_atlas, territories_labels):
    # create directory for results
    roi_stats_dir.mkdir(exist_ok=True)
    
    # set up logger
    logger_name = "HCPASL.roi_stats"
    log_out = roi_stats_dir/"roi_stats.log"
    logger = setup_logger(logger_name, log_out, "INFO")
    logger.info("Producing summary statistics in ASLT1w ROIs.")

    # get an FSL identity transform from asl2struct
    identity_name = roi_stats_dir/"fsl_identity.txt"
    rt.Registration.identity().save_fsl(str(identity_name), src=str(gm_pve), ref=str(struct_name))
    
    # set up and run oxford_asl_roi_stats command
    roi_script_name = Path(os.environ["FSLDIR"])/"bin/oxford_asl_roi_stats.py"
    cmd = [
        "fslpython", str(roi_script_name),
        "--oxasl-output", str(oxford_asl_dir),
        "--struc", str(struct_name),
        "--gm-pve", str(gm_pve), "--wm-pve", str(wm_pve),
        "--asl2struc", str(identity_name), "--std2struc", str(std2struct_name),
        "--output", str(roi_stats_dir),
        "--add-mni-atlas", str(territories_atlas), "--add-mni-atlas-labels", str(territories_labels),
        "--output-prefix", "roi_stats",
        "--gm-thresh", "0.7", "--native-pves"
    ]
    logger.info("Running oxford_asl_roi_stats.py with command:")
    logger.info(" ".join(cmd))
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    while 1:
        retcode = process.poll()
        line = process.stdout.readline().decode("utf-8")
        logger.info(line)
        if line =="" and retcode is not None:
            break
    if retcode != 0:
        logger.info(f"retcode={retcode}")
        logger.exception("Process failed.")