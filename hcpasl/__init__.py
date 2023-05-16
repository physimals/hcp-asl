"""
HCP-ASL pipeline

Copyright 2020-21 University of Nottingham
"""
try:
    from ._version import __version__, __timestamp__, __sha1__
except ImportError:
    __version__ = "unknown"
    __timestamp__ = "unknown"
    __sha1__ = "unknown"

__all__ = ["__version__", "__timestamp__", "__sha1__"]

import os 
if not os.environ.get("HCPPIPEDIR"): 
    raise RuntimeError("Environment variable HCPPIPEDIR must be set (see HCP pipeline installation)")

if not os.environ.get("FREESURFER_HOME"): 
    raise RuntimeError("Environment variable FREESURFER_HOME must be set (see FreeSurfer installation)")