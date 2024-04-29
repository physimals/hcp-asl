import io
import os
import re
import subprocess

from setuptools import find_packages, setup

PACKAGE_NAME = "hcpasl"
ROOTDIR = os.path.abspath(os.path.dirname(__file__))


def git_version():
    """Get the full and python standardized version from Git tags (if possible)"""
    try:
        # Full version includes the Git commit hash
        full_version = (
            subprocess.check_output("git describe", shell=True)
            .decode("utf-8")
            .strip(" \n")
        )

        # Python standardized version in form major.minor.patch.post<build>
        version_regex = re.compile(r"v?(\d+\.\d+\.\d+(-\d+)?).*")
        match = version_regex.match(full_version)
        if match:
            std_version = match.group(1).replace("-", ".post")
        else:
            raise RuntimeError("Failed to parse version string %s" % full_version)
        return full_version, std_version

    except Exception:
        # Any failure, return None. We may not be in a Git repo at all
        return None, None


def git_timestamp():
    """Get the last commit timestamp from Git (if possible)"""
    try:
        return (
            subprocess.check_output("git log -1 --format=%cd", shell=True)
            .decode("utf-8")
            .strip(" \n")
        )
    except Exception:
        # Any failure, return None. We may not be in a Git repo at all
        return None


def git_sha1():
    """Get the last commit SHA1 hash from Git (if possible)"""
    try:
        return (
            subprocess.check_output("git rev-parse HEAD", shell=True)
            .decode("utf-8")
            .strip(" \n")
        )
    except Exception:
        # Any failure, return None. We may not be in a Git repo at all
        return None


def update_metadata(version_str, timestamp_str, sha1_str):
    """Update the version and timestamp metadata in the module _version.py file"""
    with io.open(
        os.path.join(ROOTDIR, PACKAGE_NAME, "_version.py"), "w", encoding="utf-8"
    ) as f:
        f.write("__version__ = '%s'\n" % version_str)
        f.write("__timestamp__ = '%s'\n" % timestamp_str)
        f.write("__sha1__ = '%s'\n" % sha1_str)


def get_requirements():
    """Get a list of all entries in the requirements file"""
    with io.open(os.path.join(ROOTDIR, "requirements.txt"), encoding="utf-8") as f:
        return [l.strip() for l in f.readlines()]


def get_version():
    """Get the current version number (and update it in the module _version.py file if necessary)"""
    version, timestamp, sha1 = git_version()[1], git_timestamp(), git_sha1()

    if version is not None and timestamp is not None and sha1 is not None:
        # We got the metadata from Git - update the version file
        update_metadata(version, timestamp, sha1)
    else:
        # Could not get metadata from Git - use the version file if it already exists
        try:
            with io.open(
                os.path.join(ROOTDIR, PACKAGE_NAME, "_version.py"), encoding="utf-8"
            ) as f:
                md = f.read()
                match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", md, re.M)
                if match:
                    version = match.group(1)
                else:
                    raise ValueError("Stored version could not be parsed")
        except (IOError, ValueError):
            # No version file - create one with 'unknown' placholder data
            version = "unknown"
            update_metadata(version, "unknown", "unknown")
    return version


with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name=PACKAGE_NAME,
    version=get_version(),
    author="Flora Kennedy McConnell, Jack Toner, Thomas Kirk",
    author_email="thomas.kirk1@nottingham.ac.uk",
    description="Minimal ASL processing pipeline for the HCP Lifespan datasets",
    long_description=long_description,
    url="https://github.com/physimals/hcp-asl",
    packages=find_packages(),
    python_requires=">=3.9",
    install_requires=get_requirements(),
    entry_points={
        "console_scripts": [
            "process_hcp_asl = scripts.run_pipeline:main",
            "get_sebased_bias_asl = scripts.se_based:se_based_bias_estimation",
            "mt_estimation_asl = scripts.mt_estimation_pipeline:main",
            "results_to_mni_asl = scripts.results_to_mni:main",
        ]
    },
    scripts=[
        "scripts/VolumetoSurfaceASL.sh",
        "scripts/SurfaceSmoothASL.sh",
        "scripts/SubcorticalProcessingASL.sh",
        "scripts/PerfusionCIFTIProcessingPipelineASL.sh",
        "scripts/CreateDenseScalarASL.sh",
    ],
    package_data={
        "hcpasl": [
            "resources/empirical_banding_factors.txt",
            "resources/vascular_territories_atlas.nii.gz",
            "resources/vascular_territories_atlas_labels.txt",
            "resources/ASLQC_template.scene",
        ]
    },
)
