from setuptools import setup, find_packages

setup(
    name="martinisurf",
    version="1.0.0",
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
        "vermouth" ,
        "mdtraj",
    ],
    entry_points={
        "console_scripts": [
            "martinisurf = martinisurf.__main__:main",
        ]
    }
)
