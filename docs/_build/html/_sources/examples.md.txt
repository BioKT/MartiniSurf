# Recommended Examples

The recommended complete workflows are:

## Protein + Not explicit Linker + solvate + ionize

Path: `martinisurf/examples/protein/04_anchor_solvate_ionize`

```bash
cd martinisurf/examples/protein/04_anchor_solvate_ionize
bash run.sh
bash work_flow_gromacs.sh
```

## DNA + linker + solvate + ionize + frozen water

Path: `martinisurf/examples/dna/03_linker_solvate_ionize_freeze`

```bash
cd martinisurf/examples/dna/03_linker_solvate_ionize_freeze
bash run.sh
bash work_flow_gromacs.sh
```

Behavior notes:
- In DNA linker mode, linker-DNA is coupled with bonded terms (bond + angle) in topology (no linker-DNA pull coordinate).
- `work_flow_gromacs.sh` uses non-restrained topology in equilibration stages and restricted topology in production when available.

## DNA + linker + solvate + polarizable water

Path: `martinisurf/examples/dna/04_linker_solvate_polarizable_water`

```bash
cd martinisurf/examples/dna/04_linker_solvate_polarizable_water
bash run.sh
```

Behavior notes:
- Uses `--polarizable-water`, so the final topology includes `martini_v2.1P-dna.itp` and solvent `PW`.
- This bundled example is generated up to build + solvation. Legacy Martini 2 polarizable-water `.mdp` files may require a compatible GROMACS/MDP stack before ionization or production MD.

## Pre-CG protein+cofactor + substrate + solvate + ionize

Path: `martinisurf/examples/protein/05_pre_cg_nad_substrate`

```bash
cd martinisurf/examples/protein/05_pre_cg_nad_substrate
bash run.sh
bash work_flow_gromacs.sh
```

Behavior notes:
- In `pre_cg_complex`, low-Z balancing is enabled by default.
- Default low-Z fraction is `0.2` unless `protein.balance_low_z_fraction` is provided in `complex_config.yaml`.

## Protein build-only examples

- `martinisurf/examples/protein/01_anchor`
- `martinisurf/examples/protein/02_adsorption`
- `martinisurf/examples/protein/03_linker_surface_decoration`

## DNA build-only examples

- `martinisurf/examples/dna/01_anchor`
- `martinisurf/examples/dna/02_linker`

## Surface-focused examples

- `martinisurf/examples/surfaces/01_protein_cnt_3d`
- `martinisurf/examples/surfaces/02_protein_graphene_resizable`
- `martinisurf/examples/surfaces/03_protein_bilayer_hexagonal`

## Adsorption mode note

For adsorption setups, use the main build command with `--ads-mode` and anchor-based orientation
(`--anchor ...` or `--complex-config` anchor groups). This mode runs:
`minimization -> nvt -> npt -> production` without deposition pulls/restraints.
