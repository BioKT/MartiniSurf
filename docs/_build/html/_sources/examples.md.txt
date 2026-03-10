# Recommended Examples

The recommended complete workflows are:

## Protein + Not explicit Linker + solvate + ionize

Path: `martinisurf/examples/06_protein_anchor_solvate_ionize`

```bash
cd martinisurf/examples/06_protein_anchor_solvate_ionize
bash run.sh
bash work_flow_gromacs.sh
```

## DNA + linker + solvate + ionize + frozen water

Path: `martinisurf/examples/07_dna_linker_solvate_ionize_freeze`

```bash
cd martinisurf/examples/07_dna_linker_solvate_ionize_freeze
bash run.sh
bash work_flow_gromacs.sh
```

Behavior notes:
- In DNA linker mode, linker-DNA is coupled with bonded terms (bond + angle) in topology (no linker-DNA pull coordinate).
- `work_flow_gromacs.sh` uses non-restrained topology in equilibration stages and restricted topology in production when available.

## Pre-CG protein+cofactor + substrate + solvate + ionize

Path: `martinisurf/examples/08_protein_nad`

```bash
cd martinisurf/examples/08_protein_nad
bash run.sh
bash work_flow_gromacs.sh
```

Behavior notes:
- In `pre_cg_complex`, low-Z balancing is enabled by default.
- Default low-Z fraction is `0.2` unless `protein.balance_low_z_fraction` is provided in `complex_config.yaml`.

## Basic build-only examples

- `martinisurf/examples/01_protein_anchor`
- `martinisurf/examples/02_protein_ads`
- `martinisurf/examples/03_protein_linker`
- `martinisurf/examples/04_dna_anchor`
- `martinisurf/examples/05_dna_linker`

## Adsorption mode note

For adsorption setups, use the main build command with `--ads-mode` and anchor-based orientation
(`--anchor ...` or `--complex-config` anchor groups). This mode runs:
`minimization -> nvt -> npt -> production` without deposition pulls/restraints.
