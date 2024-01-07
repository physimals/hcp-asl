import logging
import os
from pathlib import Path
from shutil import rmtree
from string import Template

import regtricks as rt

from hcpasl.utils import get_package_data_name, get_roi_stats_script, sp_run


def create_qc_report(subject_dir, outdir, regname="MSMAll"):
    if not outdir:
        outdir = subject_dir
    else:
        outdir = subject_dir / outdir

    # Sneaky redirection of paths.
    # The template scene file is referenced to T1w/ASL/ASLQC
    # So we generate the scene file there, template in subject-specific data
    subject_id = subject_dir.stem
    scene_initial = (
        Path(subject_dir) / f"T1w/ASL/ASLQC/{subject_id}_hcp_asl_qc.wb_scene"
    )
    rel_path_outdir = os.path.relpath(outdir / "T1w/ASL", start=scene_initial.parent)
    scene_initial.parent.mkdir(parents=True, exist_ok=True)

    # Load template file
    template_name = get_package_data_name("ASLQC_template.scene")
    with open(template_name, "r") as f:
        scene_template = Template(f.read())

    # Write in subject variables
    data = {"SUBID": subject_id, "REL_PATH_TO_OUT_T1wASL": rel_path_outdir}
    with open(scene_initial, "w") as f:
        scene_file = scene_template.substitute(data)
        f.write(scene_file)

    scene_final = outdir / f"T1w/ASL/ASLQC/{subject_id}_hcp_asl_qc.wb_scene"
    if scene_final != scene_initial:
        scene_final.parent.mkdir(parents=True, exist_ok=True)
        logging.info(f"Generating wb_view QC scene at {scene_final}")
        cmd = [
            "wb_command",
            "-scene-file-relocate",
            str(scene_initial),
            str(scene_final),
        ]
        subprocess_popen(cmd)
        rmtree(scene_initial.parent)
    else:
        scene_final = scene_initial

    wb_cmd = os.environ["CARET7DIR"] + "/wb_command"
    snapdir = scene_final.parent / f"snapshots"
    snapdir.mkdir(exist_ok=True)
    for idx in range(1, 8):
        png = snapdir / f"{subject_id}_hcp_asl_qc.wb_scene{idx}.png"
        cmd = [
            wb_cmd,
            "-show-scene",
            str(scene_final),
            str(idx),
            str(png),
            "100",
            "100",
            "-use-window-size",
        ]
        sp_run(cmd)


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
    roi_script_name = get_roi_stats_script()
    cmd = [
        "fslpython",
        str(roi_script_name),
        "--oxasl-output",
        str(oxford_asl_dir / "native_space"),
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
    subprocess_popen(cmd)
