#!/usr/bin/env python3
"""
Setup script for automated glycosylation pipeline.
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read the contents of requirements.txt
requirements_file = Path(__file__).parent / "requirements.txt"
if requirements_file.exists():
    with open(requirements_file) as f:
        requirements = [line.strip() for line in f if line.strip() and not line.startswith("#")]
else:
    requirements = [
        "numpy>=1.19.0",
        "scipy>=1.5.0",
        "mdtraj>=1.9.0",
        "matplotlib>=3.3.0",
        "pandas>=1.1.0",
        "glycosylator>=0.1.0",
        "tqdm>=4.50.0",
    ]

# Read the long description
readme_file = Path(__file__).parent / "README.md"
if readme_file.exists():
    long_description = readme_file.read_text()
else:
    long_description = "Automated glycosylation pipeline for glycoprotein modeling"

setup(
    name="automated_glycosylation",
    version="1.0.0",
    author="Anacleto Silva de Souza",
    author_email="anacletosilvadesouza@usp.br",
    description="Automated glycosylation pipeline for glycoprotein structure preparation and parametrization",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT",
    packages=find_packages(),
    include_package_data=True,
    package_data={
        "automated_glycosylation": [
            "bin/*.sh",
            "scripts/**/*.py",
        ],
    },
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "glyco-prep=automated_glycosylation.cli:run_glyco_prep",
            "glyco-param=automated_glycosylation.cli:run_glyco_param",
            "glyco-orient=automated_glycosylation.cli:run_glyco_orient",
            "glyco-all=automated_glycosylation.cli:run_all_pipeline",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
    python_requires=">=3.8",
    zip_safe=False,
)
