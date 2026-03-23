# Examples

Curated set of 12 examples, now grouped by topic under `protein/`, `dna/`, and `surfaces/`. Recommended production-oriented workflows are `protein/04`, `dna/03`, and `protein/05`, plus `dna/04` as the polarizable-water build+solvate reference.

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

## Recommended Examples
1. Protein anchor + solvate + ionize  
Path: `martinisurf/examples/protein/04_anchor_solvate_ionize`
Run: `bash run.sh` then `bash work_flow_gromacs.sh`

2. DNA linker + solvate + ionize + frozen water  
Path: `martinisurf/examples/dna/03_linker_solvate_ionize_freeze`
Run: `bash run.sh` then `bash work_flow_gromacs.sh`

3. Pre-CG protein+cofactor + substrate + solvate + ionize  
Path: `martinisurf/examples/protein/05_pre_cg_nad_substrate`
Run: `bash run.sh` then `bash work_flow_gromacs.sh`

4. DNA linker + solvate + polarizable water  
Path: `martinisurf/examples/dna/04_linker_solvate_polarizable_water`
Run: `bash run.sh`
Note: this example is bundled as build + solvation because legacy Martini 2 polarizable-water MDPs may require a compatible GROMACS/MDP stack for ionization and MD.

Workflow scripts note:
- `minimization/nvt/npt/deposition` run with non-restrained topology (`system_final.top` when available).
- `production` runs with restrained topology (`system_final_res.top`) when available.

## Protein Examples
1. Protein anchor mode  
Path: `martinisurf/examples/protein/01_anchor`

2. Protein adsorption mode (`--ads-mode`)  
Path: `martinisurf/examples/protein/02_adsorption`

3. Protein linker mode + surface decoration  
Path: `martinisurf/examples/protein/03_linker_surface_decoration`

4. Protein anchor + solvate + ionize  
Path: `martinisurf/examples/protein/04_anchor_solvate_ionize`

5. Pre-CG protein + NAD + substrate  
Path: `martinisurf/examples/protein/05_pre_cg_nad_substrate`

## DNA Examples
1. DNA anchor mode  
Path: `martinisurf/examples/dna/01_anchor`

2. DNA linker mode  
Path: `martinisurf/examples/dna/02_linker`

3. DNA linker + solvate + ionize + frozen water  
Path: `martinisurf/examples/dna/03_linker_solvate_ionize_freeze`

4. DNA linker + solvate + polarizable water  
Path: `martinisurf/examples/dna/04_linker_solvate_polarizable_water`

## Surface-Focused Examples
1. Protein immobilization on a nanotube (`--surface-geometry 3d`)  
Path: `martinisurf/examples/surfaces/01_protein_cnt_3d`

2. Protein immobilization on resizable graphene  
Path: `martinisurf/examples/surfaces/02_protein_graphene_resizable`

3. Protein immobilization on a two-layer hexagonal surface  
Path: `martinisurf/examples/surfaces/03_protein_bilayer_hexagonal`

## Per-Example Layout
- `inputs/`: required input files
- `run.sh`: full command for system generation
- `README.md`: short notes

## Colab Notebooks
- `martinisurf/examples/MartiniSurf_Protein.ipynb`
- `martinisurf/examples/MartiniSurf_DNA.ipynb`
