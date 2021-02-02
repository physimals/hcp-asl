from setuptools import setup, find_packages
import os
import subprocess
import re
import io

PACKAGE_NAME = 'hcpasl'
ROOTDIR = os.path.abspath(os.path.dirname(__file__))

def git_version():
    """ Get the full and python standardized version from Git tags (if possible) """
    try:
        # Full version includes the Git commit hash
        full_version = subprocess.check_output('git describe', shell=True).decode("utf-8").strip(" \n")

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
    """ Get the last commit timestamp from Git (if possible)"""
    try:
        return subprocess.check_output('git log -1 --format=%cd', shell=True).decode("utf-8").strip(" \n")
    except Exception:
        # Any failure, return None. We may not be in a Git repo at all
        return None

def update_metadata(version_str, timestamp_str):
    """ Update the version and timestamp metadata in the module _version.py file """
    with io.open(os.path.join(ROOTDIR, PACKAGE_NAME, "_version.py"), "w", encoding='utf-8') as f:
        f.write("__version__ = '%s'\n" % version_str)
        f.write("__timestamp__ = '%s'\n" % timestamp_str)

def get_version():
    """ Get the current version number (and update it in the module _version.py file if necessary)"""
    version, timestamp = git_version()[1], git_timestamp()

    if version is not None and timestamp is not None:
        # We got the metadata from Git - update the version file
        update_metadata(version, timestamp)
    else:
        # Could not get metadata from Git - use the version file if it exists
        try:
            with io.open(os.path.join(ROOTDIR, PACKAGE_NAME, '_version.py'), encoding='utf-8') as f:
                md = f.read()
                match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", md, re.M)
                if match:
                    version = match.group(1)
                else:
                    raise ValueError("Stored version could not be parsed")
        except (IOError, ValueError):
            version = "unknown"
            update_metadata(version, "unknown")
    return version

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name=PACKAGE_NAME,
    version=get_version(),
    author='Flora Kennedy McConnell <bbzfk@exmail.nottingham.ac.uk>, Tom Kirk <thomas.kirk@eng.ox.ac.uk>, Jack Toner <bbzjt@exmail.nottingham.ac.uk>',
    author_email='jack.toner@nottingham.ac.uk',
    description='Minimal ASL processing pipeline for the HCP.',
    long_description=long_description,
    url='https://github.com/ibme-qubic/hcp-asl',
    packages=find_packages(),
    python_requires='>=3.6',
    install_requires=[
        'numpy',
        'fslpy>=3.1.0',
        'pyfab',
        'nibabel',
        'regtricks',
        'scikit-learn',
        'toblerone @ git+https://github.com/tomfrankkirk/toblerone.git',
        'gradunwarp @ git+https://github.com/Washington-University/gradunwarp.git',
        'matplotlib',
        'pandas'
    ],
    entry_points={
        'console_scripts': [
            'hcp_asl = scripts.run_pipeline:main',
            'pv_est = scripts.prepare_t1asl_space:main',
            'get_sebased_bias = scripts.se_based:se_based_bias_estimation',
            'get_updated_fabber = scripts.get_updated_fabber:main',
            'mt_estimation = scripts.mt_estimation_pipeline:main',
            'results_to_mni = scripts.results_to_mni:main',
        ]
    },
    scripts = [
        'scripts/VolumetoSurface.sh',
        'scripts/SurfaceSmooth.sh',
        'scripts/SubcorticalProcessing.sh',
        'scripts/PerfusionCIFTIProcessingPipeline.sh',
        'scripts/CreateDenseScalar.sh',
    ]
)
