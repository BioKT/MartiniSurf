from setuptools import setup, find_packages

setup(
    name="surfmartini",
    version="1.0.0",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "numpy",
        "pyvista",
        "MDAnalysis",
        "scipy",
        "vermouth" 
    ],
    entry_points={
        "console_scripts": [
            "martinisurf = surfmartini.__main__:main",
        ]
    }
)

