import os
import subprocess
from pathlib import Path

import nbformat
import regtricks as rt
from nbclient import execute
from nbparameterise import extract_parameters, parameter_values, replace_definitions

from hcpasl.utils import get_package_data_name, setup_logger


def create_qc_report(subject_dir, outdir):
    # get location of this file
    template_name = get_package_data_name("report_template.ipynb")
    template_nb = nbformat.read(template_name, as_version=4)

    # extract original parameters from the template notebook
    orig_parameters = extract_parameters(template_nb)

    # replace with parameters for this subject and execute notebook
    new_parameters = parameter_values(
        orig_parameters, subject_dir=str(subject_dir), outdir=str(outdir)
    )
    new_nb = replace_definitions(template_nb, new_parameters, execute=False)
    _ = execute(new_nb)

    # save notebook in subject's main output directory
    new_nb_name = Path(subject_dir) / outdir / "hcp_asl_report.ipynb"
    with open(new_nb_name, "w") as f:
        nbformat.write(new_nb, f)


def roi_stats(
    struct_name,
    oxford_asl_dir,
    gm_pve,
    wm_pve,
    std2struct_name,
    roi_stats_dir,
    territories_atlas,
    territories_labels,
):
    # create directory for results
    roi_stats_dir.mkdir(exist_ok=True)

    # set up logger
    logger_name = "HCPASL.roi_stats"
    log_out = roi_stats_dir / "roi_stats.log"
    logger = setup_logger(logger_name, log_out, "INFO")
    logger.info("Producing summary statistics in ASLT1w ROIs.")

    # get an FSL identity transform from asl2struct
    identity_name = roi_stats_dir / "fsl_identity.txt"
    rt.Registration.identity().save_fsl(
        str(identity_name), src=str(gm_pve), ref=str(struct_name)
    )

    # set up and run oxford_asl_roi_stats command
    roi_script_name = Path(os.environ["FSLDIR"]) / "bin/oxford_asl_roi_stats.py"
    cmd = [
        "fslpython",
        str(roi_script_name),
        "--oxasl-output",
        str(oxford_asl_dir),
        "--struc",
        str(struct_name),
        "--gm-pve",
        str(gm_pve),
        "--wm-pve",
        str(wm_pve),
        "--asl2struc",
        str(identity_name),
        "--std2struc",
        str(std2struct_name),
        "--output",
        str(roi_stats_dir),
        "--add-mni-atlas",
        str(territories_atlas),
        "--add-mni-atlas-labels",
        str(territories_labels),
        "--output-prefix",
        "roi_stats",
        "--gm-thresh",
        "0.7",
        "--native-pves",
    ]
    logger.info("Running oxford_asl_roi_stats.py with command:")
    logger.info(" ".join(cmd))
    subprocess_popen(cmd, logger)
