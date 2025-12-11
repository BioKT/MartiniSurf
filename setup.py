from setuptools import setup, find_packages

setup(
    name="martinisurf",
    version="1.0.0",
    packages=find_packages(),
    include_package_data=True,
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

