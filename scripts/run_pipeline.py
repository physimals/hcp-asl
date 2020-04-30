"""
This script performs the full minimal pre-processing ASL pipeline 
for the Human Connectome Project (HCP) ASL data.
"""

import sys.argv
from initial_bookkeeping import initial_processing
from m0_mt_correction import correct_M0
from asl_correction import hcp_asl_moco
from asl_differencing import tag_control_differencing
from asl_perfusion import run_oxford_asl
from projection import project_to_surface

# assign arguments to variables
_, subject_dir, mt_name = sys.argv
mt_name = in_args[-1] # magnetisation transfer scaling factors
