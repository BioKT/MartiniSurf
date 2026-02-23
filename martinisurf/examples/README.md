# Examples

Curated set of 7 examples. Recommended production-oriented workflows are `05`, `06`, and `07`.

## Recommended Complete Examples
1. Protein anchor + solvate + ionize  
Path: `martinisurf/examples/05_protein_anchor_solvate_ionize`

2. DNA linker + solvate + ionize + frozen water  
Path: `martinisurf/examples/06_dna_linker_solvate_ionize_freeze`

3. Pre-CG protein+cofactor + substrate + solvate + ionize  
Path: `martinisurf/examples/07_protein_nad`

Run pattern:
- `bash run.sh`
- `bash work_flow_gromacs.sh`

## Basic Build Examples
1. Protein anchor mode  
Path: `martinisurf/examples/01_protein_anchor`

2. Protein linker mode + surface decoration  
Path: `martinisurf/examples/02_protein_linker`

3. DNA anchor mode  
Path: `martinisurf/examples/03_dna_anchor`

4. DNA linker mode  
Path: `martinisurf/examples/04_dna_linker`

## Per-Example Layout
- `inputs/`: required input files
- `run.sh`: full command for system generation
- `README.md`: short notes

## Colab Notebooks
- `martinisurf/examples/colab_protein_automartini_m3.ipynb`
- `martinisurf/examples/colab_dna_automartini_m2.ipynb`
