from .initial_bookkeeping import create_dirs
from .m0_mt_correction import load_json, update_json
from pathlib import Path
import subprocess
from fsl.wrappers.flirt import applyxfm

def project_to_surface(subject_dir):
    # load subject's json
    json_dict = load_json(subject_dir)

    # perfusion calib
    pc_name = Path(json_dict['oxford_asl']) / 'struct_space/perfusion_calib.nii.gz'

    # create directory for surface results
    projection_dir = Path(json_dict['oxford_asl']) / 'SurfaceResults'
    create_dirs([projection_dir, ])
    sides = ('L', 'R')
    for side in sides:
        # surface file names
        mid_name = json_dict[f'{side}_mid']
        pial_name = json_dict[f'{side}_pial']
        white_name = json_dict[f'{side}_white']

        # save name
        savename = projection_dir / f'{side}_perfusion_calib.func.gii'
        cmd = [
            "wb_command",
            "-volume-to-surface-mapping",
            pc_name,
            mid_name,
            savename,
            "-ribbon-constrained",
            white_name,
            pial_name
        ]
        subprocess.run(cmd)