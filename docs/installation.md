# Installation

## Local installation

```bash
conda create -n martinisurf python=3.11 -y
conda activate martinisurf
pip install -r requirements.txt
pip install -e .
```

## External tools required

- `martinize2` for protein mode
- Python 2.7 for DNA mode (`martinize-dna.py`)
- GROMACS (`gmx` or `gmx_mpi`) for solvation/ionization and simulation runs

## Quick checks

```bash
which martinisurf
which martinize2
which gmx
python2.7 -V
```

## Build docs locally

```bash
pip install -r docs/requirements.txt
sphinx-build -b html docs docs/_build/html
```

Then open `docs/_build/html/index.html`.
