from pathlib import Path

from setuptools import find_packages, setup


ROOT = Path(__file__).resolve().parent
README = (ROOT / "README.md").read_text(encoding="utf-8")

setup(
    name="martinisurf",
    version="1.0.0",
    description="Toolkit for automated Martini protein and DNA surface-system setup.",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/jjimenezgar/MartiniSurf",
    project_urls={
        "Documentation": "https://github.com/jjimenezgar/MartiniSurf/tree/master/docs",
        "Source": "https://github.com/jjimenezgar/MartiniSurf",
        "Issues": "https://github.com/jjimenezgar/MartiniSurf/issues",
    },
    python_requires=">=3.9",
    keywords=[
        "martini",
        "gromacs",
        "coarse-grained",
        "molecular-dynamics",
        "biomolecular-surfaces",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
    packages=find_packages(),
    include_package_data=True,
    package_data={
        "martinisurf": [
            "system_templates/*",
            "system_itp/*",
            "martini-dna_itp/*",
            "mdp_templates/*",
            "data/*",
            "dssp/*",
            "surfaces_generator/cnt-martini/*",
            "surfaces_generator/cnt-martini/simulation/*",
            "surfaces_generator/cnt-martini/simulation/parameters/*",
            "surfaces_generator/cnt-martini/simulation/martini/*",
            "surfaces_generator/cnt-martini/simulation/martini/gromacs-parameters/*",
            "surfaces_generator/Martini3-Graphene/*",
            "surfaces_generator/Martini3-Graphene/support/*",
            "utils/dna_backmapping_files/*",
        ],
    },
    install_requires=[
        "numpy",
        "pyvista",
        "MDAnalysis",
        "scipy",
        "vermouth",
        "mdtraj",
    ],
    entry_points={
        "console_scripts": [
            "martinisurf = martinisurf.__main__:main",
        ]
    },
    zip_safe=False,
)
