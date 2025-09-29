import json
import logging
import os
import shutil
import subprocess
from importlib.resources import path as resource_path
from dataclasses import dataclass
from typing import List, Optional, Tuple
from pathlib import Path

import nibabel as nb
import numpy as np
from fsl.data import atlases
from fsl.data.image import Image
from fsl.wrappers import LOAD, fslmaths, fslroi
from fsl.wrappers.fnirt import applywarp, invwarp

from . import resources

# ASL sequence parameters
# Default ASL sequence parameters (retain as fallbacks)
ASL_SHAPE = (86, 86, 60, 86)
NTIS = 5
IBF = "tis"
TIS = [1.7, 2.2, 2.7, 3.2, 3.7]  # s
RPTS = [6, 6, 6, 10, 15]
SLICEDT = 0.059  # s
SLICEBAND = 10  # slices per band (groups). MB factor ~= NSLICES / SLICEBAND
NSLICES = 60
BOLUS = 1.5  # s
TE = 19  # ms


@dataclass
class AslParams:
    """
    Runtime-configurable ASL acquisition parameters.

    Values are resolved by priority: JSON sidecar (if present) > CLI overrides > defaults above.

    Notes:
    - TE is stored in milliseconds to match oxford_asl expectation.
    - SLICEBAND is the number of slices within a band (group), not the multiband factor.
      For multiband factor MB and NSLICES total slices, SLICEBAND = NSLICES // MB.
        - Calibration frames: this pipeline requires exactly two calibration volumes,
            taken from the end of the series. Two preceding volumes may optionally be
            discarded to preserve legacy behavior.
    """

    # derived from NIfTI header (x, y, z, t)
    asl_shape: Tuple[int, int, int, int] = ASL_SHAPE
    nslices: int = NSLICES

    # from sidecar/CLI/defaults
    te_ms: float = TE
    sliceband: int = SLICEBAND
    slicedt: float = SLICEDT
    bolus: float = BOLUS
    ibf: str = IBF
    ntis: int = NTIS
    tis: List[float] = None
    rpts: List[int] = None

    # calibration/discard policy
    calib_vols: int = 2
    tail_discard_vols: int = 2

    def validate(self):
        if self.tis is None:
            self.tis = list(TIS)
        if self.rpts is None:
            self.rpts = list(RPTS)
        if self.ntis is None:
            self.ntis = len(self.tis)
        if len(self.tis) != len(self.rpts):
            raise ValueError(
                f"Length mismatch: TIS({len(self.tis)}) vs RPTS({len(self.rpts)})"
            )
        if self.ntis != len(self.tis):
            raise ValueError(
                f"NTIS ({self.ntis}) does not match number of TIS values ({len(self.tis)})"
            )
        if self.asl_shape is None or len(self.asl_shape) != 4:
            raise ValueError("ASL shape must be 4D (x,y,z,t)")
        if self.nslices is None:
            self.nslices = int(self.asl_shape[2])
        if self.sliceband <= 0 or self.nslices % self.sliceband != 0:
            logging.warning(
                "SLICEBAND=%d does not evenly divide NSLICES=%d; slice-timing may be inaccurate",
                self.sliceband,
                self.nslices,
            )
        if self.calib_vols < 0 or self.tail_discard_vols < 0:
            raise ValueError("Calibration and discard volumes must be >= 0")
        # enforce exactly two calibration volumes
        if self.calib_vols != 2:
            raise ValueError(
                f"Exactly 2 calibration volumes are required (got {self.calib_vols})."
            )


def _parse_csv(values: Optional[str], cast):
    if values is None:
        return None
    values = values.strip()
    if values == "":
        return None
    return [cast(v) for v in values.split(",")]


def _find_sidecar(mbpcasl_path: Path) -> Optional[Path]:
    p = Path(mbpcasl_path)
    candidates = []
    if p.name.endswith(".nii.gz"):
        candidates.append(p.with_suffix("").with_suffix(".json"))
    if p.suffix == ".nii":
        candidates.append(p.with_suffix(".json"))
    candidates.append(p.parent / (p.stem.split(".")[0] + ".json"))
    for c in candidates:
        if c.exists():
            return c
    return None


def load_asl_params(
    mbpcasl_path: Path,
    *,
    ntis: Optional[int] = None,
    tis: Optional[str] = None,
    rpts: Optional[str] = None,
    bolus: Optional[float] = None,
    slicedt: Optional[float] = None,
    te_ms: Optional[float] = None,
    sliceband: Optional[int] = None,
    calib_vols: Optional[int] = None,
    tail_discard_vols: Optional[int] = None,
    ibf: Optional[str] = None,
) -> AslParams:
    """Resolve ASL acquisition parameters for this run.

    Priority: sidecar (when available) > CLI-provided values > defaults in this module.

    - TE is reported as seconds in many BIDS sidecars; converted to milliseconds here.
        - SLICEBAND is derived as NSLICES // MultibandAccelerationFactor when available.
        - ASL_SHAPE and NSLICES are always taken from the NIfTI header.
        - Exactly two calibration volumes are required; providing a different value
            for calib_vols will raise an error.
    """

    nii = nb.load(str(mbpcasl_path))
    shp = tuple(nii.shape)

    params = AslParams(asl_shape=shp, nslices=int(shp[2]))

    # sidecar JSON
    sidecar = _find_sidecar(Path(mbpcasl_path))
    sd = None
    if sidecar is not None:
        try:
            with open(sidecar, "r") as f:
                sd = json.load(f)
        except Exception as e:
            logging.warning(f"Failed to read sidecar {sidecar}: {e}")

    # te (ms)
    if te_ms is not None:
        params.te_ms = float(te_ms)
    elif sd and "EchoTime" in sd and sd["EchoTime"] is not None:
        try:
            # Sidecar EchoTime typically in seconds
            params.te_ms = float(sd["EchoTime"]) * 1000.0
        except Exception:
            logging.warning("Could not parse EchoTime from sidecar; using default")

    # SLICEBAND from sidecar if possible
    mb_factor = None
    if sd and "MultibandAccelerationFactor" in sd:
        try:
            mb_factor = int(sd["MultibandAccelerationFactor"])
        except Exception:
            mb_factor = None
    if sliceband is not None:
        params.sliceband = int(sliceband)
    elif mb_factor and params.nslices % mb_factor == 0:
        params.sliceband = params.nslices // mb_factor

    # CLI or defaults for these
    if slicedt is not None:
        params.slicedt = float(slicedt)
    if bolus is not None:
        params.bolus = float(bolus)
    if ibf is not None:
        params.ibf = str(ibf)

    # comma separated values parsing
    tis_list = _parse_csv(tis, float)
    rpts_list = _parse_csv(rpts, int)
    if tis_list is not None:
        params.tis = tis_list
    if rpts_list is not None:
        params.rpts = rpts_list
    if ntis is not None:
        params.ntis = int(ntis)

    # calibration/discard
    if calib_vols is not None:
        params.calib_vols = int(calib_vols)
        if params.calib_vols != 2:
            raise ValueError(
                f"Exactly 2 calibration volumes are required (got {params.calib_vols})."
            )
    if tail_discard_vols is not None:
        params.tail_discard_vols = int(tail_discard_vols)

    params.validate()
    return params


class ImagePath:
    """Keep track of the name and path to an image as corrections are applied to it"""

    def __init__(self, path):
        self.path = path.resolve(strict=True)
        self.stem = self.path.stem.split(".")[0]
        self.img = nb.load(self.path)

    def correct_from_image(self, dir, suffix, newimg):
        dir.mkdir(exist_ok=True, parents=True)
        if suffix:
            stem = f"{self.stem}_{suffix}"
        else:
            stem = self.stem
        path = dir / f"{stem}.nii.gz"
        data = newimg.get_fdata()
        if data.dtype.kind == "f":
            data = data.astype(np.float32)
        else:
            data = data.astype(np.int32)
        newimg = nb.nifti1.Nifti1Image(data, affine=newimg.affine, header=newimg.header)
        nb.save(newimg, path)
        return ImagePath(path)

    def correct_from_data(self, dir, suffix, newdata):
        if newdata.dtype.kind == "f":
            newdata = newdata.astype(np.float32)
        else:
            newdata = newdata.astype(np.int32)
        newimg = nb.nifti1.Nifti1Image(newdata, self.img.affine, header=self.img.header)
        return self.correct_from_image(dir, suffix, newimg)

    def save(self):
        self.img.to_filename(self.path)

    def get_fdata(self):
        return self.img.get_fdata().astype(np.float32)

    def __str__(self):
        return str(self.path)


def load_json(subject_dir):
    """
    Load json but with some error-checking to make sure it exists.
    If it doesn't exist, instruct user to run the first part of
    the pipeline first.

    Parameters
    ----------
    subject_dir : pathlib.Path
        Path to subject's base directory.

    Returns
    -------
    dict
        Dictionary containing important file paths.
    """
    json_name = subject_dir / "ASL/ASL.json"
    if json_name.exists():
        with open(json_name, "r") as infile:
            json_dict = json.load(infile)
    else:
        raise Exception(
            f"File {json_name} does not exist."
            + "Please run initial_processing() first."
        )
    return json_dict


def update_json(new_dict, old_dict):
    """
    Add key-value pairs to the dictionary.

    Add the key-value pairs in `new_dict` to `old_dict`
    and save the resulting dictionary to the json found in
    `old_dict['json_name']`.

    Parameters
    ----------
    new_dict : dict
        New key-value pairs to add to the json of important
        file names.
    old_dict : dict
        Dictionary to be updated. Also has a field containing
        the location of the json to update.
    """
    old_dict.update(new_dict)
    with open(Path(old_dict["json_name"]), "w") as fp:
        json.dump(old_dict, fp, sort_keys=True, indent=4)


def binarise(pve_name, threshold=0.5):
    """
    Binarise a PVE image given a threshold.

    The default threshold is 0.5, the threshold used by asl_reg
    to obtain a white matter segmentation from a white matter
    PVE image.
    """
    pve = Image(str(pve_name))
    seg = Image(np.where(pve.data > threshold, 1.0, 0.0), header=pve.header)
    return seg


def linear_asl_reg(calib_name, results_dir, t1_name, t1_brain_name, wm_mask):
    # run asl_reg
    aslreg_cmd = [
        "asl_reg",
        "-i",
        calib_name,
        "-o",
        results_dir,
        "-s",
        t1_name,
        "--sbet",
        t1_brain_name,
        "--tissseg",
        wm_mask,
    ]
    sp_run(aslreg_cmd)


def get_ventricular_csf_mask(fslanatdir, interpolation=3):
    """
    Get a ventricular mask in T1 space.

    Register the ventricles mask from the harvardoxford-subcortical
    2mm atlas to T1 space, using the T1_to_MNI_nonlin_coeff.nii.gz
    in the provided fsl_anat directory.

    Parameters
    ----------
    fslanatdir: pathlib.Path
        Path to an fsl_anat output directory.
    interpolation: int
        Order of interpolation to use when performing registration.
        Default is 3.
    """
    # get ventricles mask from Harv-Ox
    atlases.rescanAtlases()
    harvox = atlases.loadAtlas("harvardoxford-subcortical", resolution=2.0)
    vent_img = Image(
        harvox.data[:, :, :, 2] + harvox.data[:, :, :, 13], header=harvox.header
    )
    vent_img = fslmaths(vent_img).thr(0.1).bin().ero().run(LOAD)

    # apply inv(T1->MNI) registration
    t1_brain = (fslanatdir / "T1_restore_brain.nii.gz").resolve(strict=True)
    struct2mni = (fslanatdir / "T1_to_MNI_nonlin_coeff.nii.gz").resolve(strict=True)
    mni2struct = invwarp(str(struct2mni), str(t1_brain), LOAD)["out"]
    vent_t1_nii = applywarp(vent_img, str(t1_brain), LOAD, warp=mni2struct)["out"]
    vent_t1_nii = fslmaths(vent_t1_nii).thr(0.9).bin().run(LOAD)
    return vent_t1_nii


def compute_split_indices(
    total_vols: int, calib_vols: int = 2, tail_discard_vols: int = 2
):
    """Compute indices for splitting ASL into data and calibrations.

    Assumes the last two volumes are calibrations and two immediately preceding
    volumes may be discarded (legacy HCP behavior). Exactly two calibration
    volumes are required.

    Returns (data_start, data_len, calib_indices_list)
    """
    if total_vols <= 0:
        raise ValueError("ASL sequence has no timepoints")
    if calib_vols < 0 or tail_discard_vols < 0:
        raise ValueError("calib_vols and tail_discard_vols must be non-negative")
    if calib_vols != 2:
        raise ValueError(
            f"Exactly 2 calibration volumes are required (got {calib_vols})."
        )
    trail = calib_vols + tail_discard_vols
    if total_vols <= trail:
        raise ValueError(
            f"Not enough volumes to split: total={total_vols}, trailing_non_data={trail}"
        )
    data_start = 0
    data_len = total_vols - trail
    calib_start = total_vols - calib_vols
    calib_idxs = list(range(calib_start, total_vols))
    return data_start, data_len, calib_idxs


def split_asl(
    asl, tis_name, calib0_name, calib1_name, *, params: Optional[AslParams] = None
):
    """
    Split ASL sequence into label-control series and calibration images.

    The split strategy is determined from the ASL image length and AslParams
    (calib_vols, tail_discard_vols). Exactly two calibration volumes are required;
    additional or fewer calibration volumes are not supported.
    """
    asl_img = nb.load(str(asl))
    total_vols = asl_img.shape[3] if asl_img.ndim == 4 else 1
    if params is None:
        # default behavior: 2 calib, 2 discarded
        calib_vols = 2
        tail_discard_vols = 2
    else:
        calib_vols = params.calib_vols
        tail_discard_vols = params.tail_discard_vols

    _, data_len, calib_idxs = compute_split_indices(
        total_vols, calib_vols=calib_vols, tail_discard_vols=tail_discard_vols
    )
    # cata
    fslroi(str(asl), str(tis_name), 0, data_len)

    # calibrations (require exactly two)
    if len(calib_idxs) != 2:
        raise ValueError(
            f"Expected exactly two calibration volumes at end of series (got {len(calib_idxs)})."
        )
    fslroi(str(asl), str(calib0_name), calib_idxs[0], 1)
    fslroi(str(asl), str(calib1_name), calib_idxs[1], 1)


def setup_logger(file_path):
    """
    Convenience function which returns a logger object with a
    specified name and level of reporting.

    Parameters
    ----------
    logger_name : str
        Name of the logger
    out_name : str
        Desired path for the logfile
    level : str (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        Desired level of reporting for the log. See Python's logging
        documentation for more info.
    verbose : bool, default=False
        If True, information will also be sent to the terminal via a
        StreamHandler as well as to a log file.
    mode : str, default="w"
        The mode of operation for the FileHandler. The default mode,
        "w", overwrites a logfile of the same name if it exists.
    """

    # set up logger's base reporting level and formatting
    logger = logging.getLogger()
    logger.setLevel("INFO")
    formatter = logging.Formatter(
        "%(asctime)s %(levelname)s %(module)s/%(funcName)s: %(message)s"
    )

    # set up FileHandler and StreamHanlder
    handlers = [logging.FileHandler(file_path, mode="w")]
    handlers.append(logging.StreamHandler())

    # add formatting to handlers and add to logger
    for handler in handlers:
        handler.setFormatter(formatter)
        logger.addHandler(handler)


def get_package_data_name(name):
    """
    This function returns the filepath to a given file which are
    included as part of the pipeline.
    """
    p = resource_path(resources, name)
    with p as filename:
        name = filename
    return name


def copy_oxford_asl_inputs(input_dict, output_dir):
    """
    Take a dictionary of inputs to oxford_asl and copy them into one directory.
    """
    # create new output directory
    output_dir.mkdir(exist_ok=True, parents=True)

    # create dict of new filenames
    new_dict = {key: output_dir / val.name for key, val in input_dict.items()}

    # copy files from original locations in input_dict to new location in new_dict
    for old_loc, new_loc in zip(input_dict.values(), new_dict.values()):
        shutil.copy2(old_loc, new_loc)


def make_motion_fov_mask(mc_transform, src, ref, cores=1):
    """Generate a FoV array in the reference of the motion transform
    Invert the transform and map the FoV into the source space
    Apply logical all across time to get valid FoV voxels
    Save mask in source space"""

    fov = np.ones(ref.size)
    fov_motion = mc_transform.apply_to_array(fov, ref, src, order=1, cores=cores)
    fov_valid = (fov_motion > 0.9).all(-1)
    return src.make_nifti(fov_valid)


def sp_run(cmd, **kwargs):
    logging.info(cmd)
    env = {**os.environ, **kwargs.pop("env", {})}
    result = subprocess.run(cmd, capture_output=True, text=True, env=env, **kwargs)
    if result.returncode == 0:
        logging.info(result.stdout)
    else:
        logging.error(f"Subprocess {cmd} failed with exit code {result.returncode}.")
        logging.error(result.stderr)
        exit(-1)


def get_roi_stats_script():
    """FSL script can have two different names (with or without .py extension)"""
    roi_script_name = Path(os.environ["FSLDIR"]) / "bin/oxford_asl_roi_stats"
    if not roi_script_name.exists():
        roi_script_name = roi_script_name.with_suffix(".py")
    if not roi_script_name.exists():
        raise RuntimeError("Cannot find oxford_asl_roi_stats within $FSLDIR")
    return roi_script_name
