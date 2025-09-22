import logging
import os.path as op

import numpy as np
import regtricks as rt

from .utils import sp_run


def register_asl2struct(src, struct, fsdir, reg_dir):
    """
    Generate the linear transformation between ASL-space and T1w-space
    using FS bbregister. Note that struct is required only for saving
    the output in the right convention, it is not actually used by
    bbregister.

    Args:
        src: path to volume in ASL voxel grid
        struct: path to T1w image (eg T1w_acdc_restore.nii.gz)
        fsdir: path to subject's FreeSurfer output directory
        reg_dir: path to registration directory, for output

    Returns:
        n/a, file 'asl2struct.mat' will be saved in reg_dir
    """

    logging.info(f"Movable volume: {src}")
    logging.info(f"T1w structural image: {struct}")
    logging.info(f"FreeSurfer output directory: {fsdir}")
    logging.info(f"Output directory: {reg_dir}")

    # We need to do some hacky stuff to get bbregister to work...
    # Split the path to the FS directory into a fake $SUBJECTS_DIR
    # and subject_id. We temporarily set the environment variable
    # before making the call, and then revert back afterwards
    new_sd, sid = op.split(fsdir)
    orig_mgz = op.join(fsdir, "mri", "orig.mgz")

    # Save the output in fsl format, by default
    # this targets the orig.mgz, NOT THE T1 IMAGE ITSELF!
    logging.info(f"Running bbregister in {reg_dir}; setting $SUBJECTS_DIR to {new_sd}")
    omat_path = op.join(reg_dir, "asl2struct.mat")
    cmd = f"$FREESURFER_HOME/bin/bbregister --s {sid} --mov {src} --t2 "
    cmd += f"--reg asl2orig_mgz_initial_bbr.dat --fslmat {omat_path} --init-fsl"
    fslog_name = op.join(reg_dir, "asl2orig_mgz_initial_bbr.dat.log")
    logging.info(f"FreeSurfer's bbregister log: {fslog_name}")
    sp_run(cmd, shell=True, env={"SUBJECTS_DIR": new_sd}, cwd=reg_dir)

    # log final .dat transform
    with open(omat_path, "r") as f:
        lines = f.readlines()
        for line in lines:
            logging.info(line)

    # log minimum registration cost
    mincost = np.loadtxt(op.join(reg_dir, "asl2orig_mgz_initial_bbr.dat.mincost"))
    logging.info(f"bbregister's mincost: {mincost[0]:4f}")

    # convert .dat to .mat
    try:
        asl2orig_fsl = rt.Registration.from_flirt(omat_path, src, orig_mgz)
    except RuntimeError:
        # final row != [0 0 0 1], round to 5 d.p. and try again
        logging.warning("FSL .mat file has an invalid format. Rounding to 5 d.p.")
        arr = np.loadtxt(omat_path)
        np.savetxt(omat_path, arr, fmt="%.5f")
        asl2orig_fsl = rt.Registration.from_flirt(omat_path, src, orig_mgz)

    # Return to original working directory, and flip the FSL matrix to target
    # asl -> T1, not orig.mgz. Save output.
    logging.info("Converting .mat to target T1w.nii.gz rather than orig.mgz")
    asl2struct_fsl = asl2orig_fsl.to_flirt(src, struct)
    np.savetxt(op.join(reg_dir, "asl2struct.mat"), asl2struct_fsl)
