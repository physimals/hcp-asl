import argparse
from pathlib import Path
from setup_mt_estimation import setup_mtestimation
from estimate_MT import estimate_mt
import multiprocessing as mp
from functools import partial

TR = 8
ROIS = {
    "wm": ["wm"],
    "gm": ["gm"],
    "csf": ["csf"],
    "combined": ["combined"],
    "all": ["wm", "gm", "csf", "combined"]
}

def main():
    # argument handling
    parser = argparse.ArgumentParser(
        description="Run the MT estimation pipeline."
    )
    parser.add_argument("--subjectlist",
        help="A .txt file of subject names from whom we "
            +"wish to estimate the MT scaling factors.",
        required=True
    )
    parser.add_argument("--studydir",
        help="Path to the study's base directory.",
        required=True
    )
    parser.add_argument("--method",
        help="Method of bias correction to use. Default is "
            +"'calib'.",
        default='calib',
        choices=('calib', 't1', 'sebased')
    )
    parser.add_argument("-c", "--cores",
        help="Number of cores to use. Default is the number "
            +f"of cores your machine has ({mp.cpu_count()}).",
        default=mp.cpu_count(),
        type=int,
        choices=range(1, mp.cpu_count()+1)
    )
    parser.add_argument("--roi",
        help="Tissue in which to estimate the MT scaling factors.",
        default="combined",
        choices=("combined", "wm", "gm", "csf", "all")
    )
    parser.add_argument("-v", "--verbose",
        help="Print some useful statements.",
        action="store_true"
    )

    # parse
    studydir = Path(args.studydir).resolve(strict=True)
    subjects = Path(args.subjects).resolve(strict=True)
    subjects = np.loadtxt(subjects, dtype=str)
    subjects = [studydir / subid for subid in subjects]
    rois = ROIS[args.roi]
    if args.verbose:
        print(f"Your study directory is {studydir}.")
        print(f"You are processing {len(subjects)} subjects.")
        print(f"Rois to be used for the estimation: {rois}.")

    # do setup
    setup_call = partial(setup_mtestimation, rois)
    with mp.Pool(args.cores) as pool:
        results = pool.map(setup_call, subjects)
    for result in results:
        print(result)

    # do estimation
    estimate_mt(subject_dirs, rois, TR, 'separate')