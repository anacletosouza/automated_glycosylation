from setuptools import setup, find_packages

setup(
    name="automated_glycosylation",
    version="1.0.0",
    author="Anacleto Silva de Souza",
    author_email="anacletosilvadesouza@usp.br",
    description="Automated glycosylation pipeline for glycoproteins",
    license="MIT",
    packages=find_packages(),
    include_package_data=True,
    package_data={
        "automated_glycosylation": [
            "scripts/**/*.py",
            "scripts/1_glycosylation_preparation/*.py",
            "scripts/2_parametrization_scripts/*.py",
            "scripts/3_carbohydrate_orientation/*.py",
        ],
    },
    entry_points={
        "console_scripts": [
            "glyco-prep=automated_glycosylation.cli:run_glyco_prep",
            "glyco-param=automated_glycosylation.cli:run_glyco_param",
            "glyco-orient=automated_glycosylation.cli:run_glyco_orient",
            "glyco-all=automated_glycosylation.cli:run_all_pipeline",
        ],
    },
    install_requires=[
        "numpy>=1.19.0",
        # Adicione outras dependências necessárias
    ],
    python_requires=">=3.7",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
