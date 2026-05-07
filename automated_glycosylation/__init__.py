#!/usr/bin/env python3
"""
Automated Glycosylation Pipeline

A comprehensive pipeline for automated glycosylation of proteins,
including glycan attachment, parametrization, and orientation optimization.
"""

from .cli import auto_glyco

__version__ = "1.0.0"
__author__ = "Anacleto Silva de Souza"
__email__ = "anacletosilvadesouza@usp.br"

__all__ = [
    "auto_glyco",
    "__version__",
    "__author__",
    "__email__",
]
