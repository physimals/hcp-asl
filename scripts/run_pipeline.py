"""
This script performs the full minimal pre-processing ASL pipeline 
for the Human Connectome Project (HCP) ASL data.

This currently requires that the script is called followed by 
the directories of the subjects of interest and finally the 
name of the MT correction scaling factors image.
"""

import sys.argv
from hcp_asl.initial_bookkeeping import initial_processing
from hcp_asl.m0_mt_correction import correct_M0
from hcp_asl.asl_correction import hcp_asl_moco
from hcp_asl.asl_differencing import tag_control_differencing
from hcp_asl.asl_perfusion import run_oxford_asl
from hcp_asl.projection import project_to_surface
from pathlib import Path
import subprocess

# assign arguments to variables
_, subject_dir, mt_name = sys.argv
mt_name = in_args[-1] # magnetisation transfer scaling factors

def main(subject_dir, mt_factors):
    """
    Run pipeline for individual subject specified by 
    `subject_dir`.
    """
    subject_dir = Path(subject_dir)
    mt_factors = Path(mt_factors)
    initial_processing(subject_dir)
    correct_M0(subject_dir, mt_factors)
    hcp_asl_moco(subject_dir, mt_factors)
    dist_corr_call = [
        "distcorr_warps",
        subject_dir.parent,
        subject_dir.stem
    ]
    subprocess.run(dist_corr_call, check=True)
    tag_control_differencing(subject_dir)
    run_oxford_asl(subject_dir)
    project_to_surface(subject_dir)

if __name__ == '__main__':
    # assign arguments to variables
    in_args = sys.argv[1:]
    mt_name = in_args[-1]
    subject_list = in_args[:-1]
    for subject_dir in subject_list:
        main(subject_dir)