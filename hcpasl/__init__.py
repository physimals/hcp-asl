"""
HCP-ASL pipeline

Copyright 2020-25 University of Nottingham, Quantified Imaging
"""

try:
    from ._version import __sha1__, __timestamp__, __version__
except ImportError:
    __version__ = "unknown"
    __timestamp__ = "unknown"
    __sha1__ = "unknown"

__all__ = ["__version__", "__timestamp__", "__sha1__"]
