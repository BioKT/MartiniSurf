<p align="center">
  <img src="logo_martini.png" alt="MartiniSurf Logo">
</p>

<h1 align="center">MartiniSurf</h1>

<p align="center">
Toolkit for automated Martini protein/DNA surface-system setup, including linker-aware orientation and pull/bonded coupling generation.
</p>

<p align="center">
  <a href="https://github.com/BioKT/MartiniSurf/actions/workflows/python-ci.yml">
    <img src="https://github.com/BioKT/MartiniSurf/actions/workflows/python-ci.yml/badge.svg" alt="CI">
  </a>
  <a href="https://biokt.github.io/MartiniSurf/">
    <img src="https://img.shields.io/badge/docs-GitHub%20Pages-2ea44f?logo=github" alt="Docs">
  </a>
  <img src="https://img.shields.io/badge/python-3.9%20|%203.10%20|%203.11-blue.svg" alt="Python versions">
  <a href="https://colab.research.google.com/github/BioKT/MartiniSurf/blob/master/martinisurf/examples/MartiniSurf_Protein.ipynb">
    <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open MartiniSurf_Protein">
  </a>
  <a href="https://colab.research.google.com/github/BioKT/MartiniSurf/blob/master/martinisurf/examples/MartiniSurf_DNA.ipynb">
    <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open MartiniSurf_DNA">
  </a>
</p>

## Overview
MartiniSurf builds complete GROMACS-ready simulation folders for:
- Protein-surface systems (Martini 3 workflow)
- DNA-surface systems (Martini2 DNA workflow)
- Linker-mediated immobilization workflows

## Documentation
- Published docs (GitHub Pages): https://biokt.github.io/MartiniSurf/
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


## Quick Start (Recommended Examples)
Examples are grouped by topic under `martinisurf/examples/protein`, `martinisurf/examples/dna`, and `martinisurf/examples/surfaces`.
The main ready-to-use workflows are:

- `martinisurf/examples/protein/04_anchor_solvate_ionize`
- `martinisurf/examples/dna/03_linker_solvate_ionize_freeze`
- `martinisurf/examples/protein/05_pre_cg_nad_substrate`

Typical run pattern for the full workflow examples (`protein/04`, `dna/03`, `protein/05`):
```bash
cd martinisurf/examples/protein/04_anchor_solvate_ionize
bash run.sh
bash work_flow_gromacs.sh
```

DNA workflow note:
- `dna/03` runs `minimization -> nvt -> deposition (NPT) -> production (NVT)`.
- The pressure-coupled equilibration is handled in `deposition`; there is no separate DNA `npt.mdp` stage in this example.

## Google Colab Notebooks
- Protein workflow + optional linker generation with AutoMartini M3:
  - Notebook: `martinisurf/examples/MartiniSurf_Protein.ipynb`
  - Open in Colab: https://colab.research.google.com/github/BioKT/MartiniSurf/blob/master/martinisurf/examples/MartiniSurf_Protein.ipynb
- DNA workflow + optional linker generation with auto_martini M2:
  - Notebook: `martinisurf/examples/MartiniSurf_DNA.ipynb`
  - Open in Colab: https://colab.research.google.com/github/BioKT/MartiniSurf/blob/master/martinisurf/examples/MartiniSurf_DNA.ipynb

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
