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
