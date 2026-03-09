# Recommended Examples

The recommended complete workflows are:

## Protein + Not explicit Linker + solvate + ionize

Path: `martinisurf/examples/05_protein_anchor_solvate_ionize`

```bash
cd martinisurf/examples/05_protein_anchor_solvate_ionize
bash run.sh
bash work_flow_gromacs.sh
```

## DNA + linker + solvate + ionize + frozen water

Path: `martinisurf/examples/06_dna_linker_solvate_ionize_freeze`

```bash
cd martinisurf/examples/06_dna_linker_solvate_ionize_freeze
bash run.sh
bash work_flow_gromacs.sh
```

## Pre-CG protein+cofactor + substrate + solvate + ionize

Path: `martinisurf/examples/07_protein_nad`

```bash
cd martinisurf/examples/07_protein_nad
bash run.sh
bash work_flow_gromacs.sh
```

## Basic build-only examples

- `martinisurf/examples/01_protein_anchor`
- `martinisurf/examples/02_protein_linker`
- `martinisurf/examples/03_dna_anchor`
- `martinisurf/examples/04_dna_linker`

## Adsorption mode note

For adsorption setups, use the main build command with `--ads-mode` and anchor-based orientation
(`--anchor ...` or `--complex-config` anchor groups). This mode runs:
`minimization -> nvt -> npt -> production` without deposition pulls/restraints.
