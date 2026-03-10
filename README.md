<p align="center">
  <img src="logo_martini.png" alt="MartiniSurf Logo">
</p>

<h1 align="center">MartiniSurf</h1>

<p align="center">
Toolkit for automated Martini protein/DNA surface-system setup, including linker-aware orientation and pull/bonded coupling generation.
</p>

<p align="center">
  <a href="https://github.com/jjimenezgar/MartiniSurf/actions/workflows/python-ci.yml">
    <img src="https://github.com/jjimenezgar/MartiniSurf/actions/workflows/python-ci.yml/badge.svg" alt="CI">
  </a>
  <a href="https://codecov.io/gh/jjimenezgar/MartiniSurf">
    <img src="https://codecov.io/gh/jjimenezgar/MartiniSurf/branch/master/graph/badge.svg" alt="Coverage">
  </a>
  <img src="https://img.shields.io/badge/python-3.9%20|%203.10%20|%203.11-blue.svg" alt="Python versions">
  <a href="https://colab.research.google.com/github/jjimenezgar/MartiniSurf/blob/master/martinisurf/examples/MartiniSurf_Protein.ipynb">
    <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open MartiniSurf_Protein">
  </a>
  <a href="https://colab.research.google.com/github/jjimenezgar/MartiniSurf/blob/master/martinisurf/examples/MartiniSurf_DNA.ipynb">
    <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open MartiniSurf_DNA">
  </a>
</p>

## Overview
MartiniSurf builds complete GROMACS-ready simulation folders for:
- Protein-surface systems (Martini 3 workflow)
- DNA-surface systems (Martini2 DNA workflow)
- Linker-mediated immobilization workflows

## Documentation
- Beginner-friendly complete guide (including full flag-by-flag reference): `docs/USER_GUIDE.md`

Main capabilities:
- Coarse-graining via `martinize2` (protein) or `martinize-dna.py` (DNA)
- Surface generation or reuse of provided surfaces
- Not explicit linker (anchor) or linker-based orientation
- Automatic topology assembly

## Installation
```bash
conda create -n martinisurf python=3.11 -y
conda activate martinisurf
pip install -r requirements.txt
pip install -e .
```

## External Tools
MartiniSurf expects the following tools in your environment:
- `martinize2` for protein mode
- Python 2.7 for DNA mode (`martinize-dna.py`)
- GROMACS for running the generated workflows


## Quick Start (Recommended Complete Examples)
These are the most complete examples (system build + solvation/ionization + simulation workflow):

- `martinisurf/examples/06_protein_anchor_solvate_ionize`
- `martinisurf/examples/07_dna_linker_solvate_ionize_freeze`
- `martinisurf/examples/08_protein_nad`

Run any of them with:
```bash
cd martinisurf/examples/06_protein_anchor_solvate_ionize
bash run.sh
bash work_flow_gromacs.sh
```

## Google Colab Notebooks
- Protein workflow + optional linker generation with AutoMartini M3:
  - Notebook: `martinisurf/examples/MartiniSurf_Protein.ipynb`
  - Open in Colab: https://colab.research.google.com/github/jjimenezgar/MartiniSurf/blob/master/martinisurf/examples/MartiniSurf_Protein.ipynb
- DNA workflow + optional linker generation with auto_martini M2:
  - Notebook: `martinisurf/examples/MartiniSurf_DNA.ipynb`
  - Open in Colab: https://colab.research.google.com/github/jjimenezgar/MartiniSurf/blob/master/martinisurf/examples/MartiniSurf_DNA.ipynb

## Linker Mode Notes
- `--anchor` and `--linker-group` accept either legacy global residue ids or chain-based residue ids from the input PDB:
  - Legacy: `--anchor 1 8 10 11`
  - Chain-based: `--anchor B 8 10 11`
  - Chain-based groups are converted internally to global residue ids in input order, so the first group still becomes `Anchor_1`, the second `Anchor_2`, and so on.
- Chain-based syntax is available for `--pdb` workflows.
- In `--complex-config`, chain-based `protein.anchor_groups` are also supported when `protein.reference_pdb` points to the source PDB used to build the pre-CG complex.
- The linker topology file must exist next to linker GRO with matching basename:
  - Example: `linker.gro` -> `linker.itp`
- You can reverse linker orientation with:
  - `--invert-linker`
- For multiple linker instances, pull/index groups are generated per linker automatically.
- If not provided manually, linker distances are estimated from Martini bead-size sigma rules.

## CLI Help
Use:
```bash
martinisurf -h
```
The help output is grouped by blocks:
- Input and molecule
- Martinization controls
- Surface controls
- Classical anchor mode
- Linker mode
- Output

## Output Structure
By default, MartiniSurf writes:
```text
Simulation_Files/
  0_topology/
    system.top
    system_res.top
    index.ndx
    system_itp/
  1_mdp/
  2_system/
```

## Testing
Run test suite:
```bash
pytest -q
```

## Third-Party Licensing Notice
MartiniSurf interfaces with external scientific tools and libraries that keep their original licenses.
You are responsible for complying with licenses of dependencies and external binaries in your environment.
