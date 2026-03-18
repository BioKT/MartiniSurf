# Examples

Curated set of 11 examples. Recommended production-oriented workflows are `06`, `07`, and `08`.

For examples driven from `--pdb`, the recommended syntax is now chain-based:
- `--anchor A 8 10 11`
- `--linker-group D 8 10 11`

MartiniSurf resolves those chain-local residues to the internal global residue ids automatically.

For `complex_config.yaml`, the equivalent chain-based syntax is available through `protein.anchor_groups` when `protein.reference_pdb` is provided.
In `pre_cg_complex`, low-Z balancing is enabled by default with `balance_low_z_fraction=0.2` unless overridden.

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
Path: `martinisurf/examples/06_protein_anchor_solvate_ionize`

2. DNA linker + solvate + ionize + frozen water  
Path: `martinisurf/examples/07_dna_linker_solvate_ionize_freeze`

3. Pre-CG protein+cofactor + substrate + solvate + ionize  
Path: `martinisurf/examples/08_protein_nad`

Run pattern:
- `bash run.sh`
- `bash work_flow_gromacs.sh`

Workflow scripts note:
- `minimization/nvt/npt/deposition` run with non-restrained topology (`system_final.top` when available).
- `production` runs with restrained topology (`system_final_res.top`) when available.

## Basic Build Examples
1. Protein anchor mode  
Path: `martinisurf/examples/01_protein_anchor`

2. Protein adsorption mode (`--ads-mode`)  
Path: `martinisurf/examples/02_protein_ads`

3. Protein linker mode + surface decoration  
Path: `martinisurf/examples/03_protein_linker`

4. DNA anchor mode  
Path: `martinisurf/examples/04_dna_anchor`

5. DNA linker mode  
Path: `martinisurf/examples/05_dna_linker`

6. Protein immobilization on a nanotube (`--surface-geometry 3d`)  
Path: `martinisurf/examples/09_protein_cnt_3d`

7. Protein immobilization on resizable graphene  
Path: `martinisurf/examples/10_protein_graphene_resizable`

8. Protein immobilization on a two-layer hexagonal surface  
Path: `martinisurf/examples/11_protein_bilayer_hexagonal`

## Per-Example Layout
- `inputs/`: required input files
- `run.sh`: full command for system generation
- `README.md`: short notes

## Colab Notebooks
- `martinisurf/examples/MartiniSurf_Protein.ipynb`
- `martinisurf/examples/MartiniSurf_DNA.ipynb`
