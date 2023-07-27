import os
from pathlib import Path
import logging
from string import Template


import nbformat
import regtricks as rt
from nbclient import execute
from nbparameterise import extract_parameters, parameter_values, replace_definitions

from hcpasl.utils import get_package_data_name, subprocess_popen


def create_qc_report(subject_dir, outdir):
    if not outdir:
        outdir = subject_dir
    else:
        outdir = subject_dir / outdir

    # Sneaky redirection of paths.
    # HCPASL outputs need to be referenced to the output directory
    # The scene file is inserted into the subjects T1w preprocessed directory
    # (not necessarily the same as the output) so it can pick up all the surfaces
    # etc automatically. Then the pipeline outputs are templated in using
    # relative paths. Finally we re-base the whole scene using wb_command.
    subject_id = subject_dir.stem
    scene_initial = Path(subject_dir) / f"T1w/ASL/{subject_id}_hcp_asl_qc.scene"
    rel_path_outdir = os.path.relpath(outdir / "T1w/ASL", start=scene_initial.parent)
    scene_initial.parent.mkdir(parents=True, exist_ok=True)

    # Load template file
    template_name = get_package_data_name("ASLQC_template.scene")
    with open(template_name, "r") as f:
        scene_template = Template(f.read())

    # Write in subject variables, save in T1w/ASL
    data = {"SUBID": subject_id, "REL_PATH_TO_OUT_T1wASL": rel_path_outdir}
    with open(scene_initial, "w") as f:
        scene_file = scene_template.substitute(data)
        f.write(scene_file)

    # Rebase scene to MNINonLinear/ASL/ASLQC, remove the one at T1w/ASL
    scene_final = (
        Path(subject_dir)
        / outdir
        / f"MNINonLinear/ASL/ASLQC/{subject_id}_hcp_asl_qc.scene"
    )
    scene_final.parent.mkdir(parents=True, exist_ok=True)
    logging.info(f"Generating wb_view QC scene at {scene_final}")
    cmd = ["wb_command", "-scene-file-relocate", scene_initial, scene_final]
    subprocess_popen(cmd)
    os.remove(scene_initial)

    wb_cmd = os.environ["CARET7DIR"] + "/wb_command"
    for idx in range(1, 8):
        png = scene_final.parent / f"{subject_id}_hcp_asl_qc_scene_{idx}.png"
        cmd = [
            wb_cmd,
            "-show-scene",
            str(scene_final),
            str(idx),
            png,
            "1400",
            "600",
        ]
        subprocess_popen(cmd)

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
    new_nb_name = (
        Path(subject_dir)
        / outdir
        / f"MNINonLinear/ASL/ASLQC/{subject_id}_hcp_asl_report.ipynb"
    )
    logging.info(f"Generating jupyter notebook at {new_nb_name}")
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
    roi_stats_dir.mkdir(exist_ok=True, parents=True)

    logging.info("Producing summary statistics in ASLT1w ROIs.")

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
    logging.info("Running oxford_asl_roi_stats.py with command:")
    logging.info(" ".join(cmd))
    subprocess_popen(cmd)
