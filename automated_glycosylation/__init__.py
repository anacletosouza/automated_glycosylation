"""
Automated Glycosylation Pipeline for Glycoproteins

A comprehensive pipeline for automated glycosylation of proteins,
including parametrization and carbohydrate orientation optimization.
"""

__version__ = "1.0.0"
__author__ = "Anacleto Silva de Souza"
__email__ = "anacletosilvadesouza@usp.br"

from .cli import (
    run_glyco_prep,
    run_glyco_param,
    run_glyco_orient,
    run_all_pipeline
)

__all__ = [
    "run_glyco_prep",
    "run_glyco_param",
    "run_glyco_orient",
    "run_all_pipeline",
    "__version__",
    "__author__",
    "__email__"
]
