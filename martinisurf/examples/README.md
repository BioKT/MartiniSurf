# Examples

Curated set of 7 examples. Recommended production-oriented workflows are `05`, `06`, and `07`.

For examples driven from `--pdb`, the recommended syntax is now chain-based:
- `--anchor A 8 10 11`
- `--linker-group D 8 10 11`

MartiniSurf resolves those chain-local residues to the internal global residue ids automatically.

For `complex_config.yaml`, the equivalent chain-based syntax is available through `protein.anchor_groups` when `protein.reference_pdb` is provided.

## Adsorption Mode (`--ads-mode`)

Use `--ads-mode` when you want the same orientation behavior as classical anchor mode, but without generating anchor pulling/restraint topology.

Behavior:
- Uses anchor-based orientation (`--anchor ...` or `--complex-config` with `protein.anchor_groups`/`protein.orient_by_residues`).
- Skips anchor pull definitions in MDP files.
- Skips anchor restrained topology usage in final top files.
- Keeps standard simulation stages as:
  `minimization -> nvt -> npt -> production`
- Does not generate/use `deposition` stage in this mode.

Constraints:
- `--ads-mode` is incompatible with linker mode (`--linker` / `--use-linker`).

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
- `martinisurf/examples/MartiniSurf_Protein.ipynb`
- `martinisurf/examples/MartiniSurf_DNA.ipynb`
