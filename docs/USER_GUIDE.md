# MartiniSurf Quick Guide for New Users

This guide is for users with no prior Martini or coarse-grained simulation experience.

Goal:
- Know which command to run.
- Understand what each flag does.
- Avoid common setup errors.

## 1) 30-second rules

- Choose only ONE orientation mode:
  - `--anchor ...` (classic anchor mode), or
  - `--linker ... --linker-group ...` (linker mode).
- If you do NOT pass `--surface`, you must pass `--lx` and `--ly`.
- `--ionize` always requires `--solvate`.
- `--freeze-water-fraction` works only with `--dna` and `--solvate`.
- If you already have a CG system (protein + cofactor), use `--complex-config`.

## 2) Minimum installation

```bash
pip install -r requirements.txt
pip install -e .
```

External tools:
- Protein mode: `martinize2`
- DNA mode: Python 2.7
- Solvation/ionization: `gmx` (GROMACS)

Quick checks:

```bash
which martinisurf
which martinize2
which gmx
python2.7 -V
```

## 3) Pick your workflow

- Atomistic protein PDB to surface:
use protein mode + `--anchor`.
- Protein with linker to surface:
use protein mode + `--linker` + `--linker-group`.
- Atomistic DNA PDB to surface:
use `--dna` + `--anchor` or `--linker`.
- Already CG protein+cofactor complex:
use `--complex-config`.

## 4) Copy-and-run commands

### A) Protein + anchors (best first run)

```bash
martinisurf \
  --pdb 1RJW \
  --moltype Protein \
  --surface-mode 4-1 \
  --lx 15 --ly 15 --dx 0.47 \
  --surface-bead P4 \
  --anchor 1 8 10 11 \
  --anchor 2 1025 1027 1028 \
  --merge A,B,C,D
```

### B) Protein + linker

```bash
martinisurf \
  --pdb 1RJW \
  --moltype Protein \
  --surface-mode 4-1 \
  --lx 15 --ly 15 --dx 0.47 \
  --surface-bead P4 \
  --linker inputs/EPOXY.gro \
  --linker-group 1 8 10 11 \
  --linker-group 2 1025 1027 1028 \
  --merge A,B,C,D
```

Note: if you use `EPOXY.gro`, `EPOXY.itp` must exist next to it.

### C) DNA + anchors

```bash
martinisurf \
  --dna \
  --pdb 4C64.pdb \
  --surface inputs/surface.gro \
  --anchor 1 1 \
  --dist 10 \
  --merge A,B
```

### D) DNA + linker + solvation + ions

```bash
martinisurf \
  --dna \
  --pdb 4C64.pdb \
  --surface-mode 4-1 \
  --lx 10 --ly 10 --dx 0.27 \
  --surface-bead C1 \
  --linker inputs/ALK.gro \
  --linker-group 1 1 \
  --solvate \
  --ionize \
  --salt-conc 0.15 \
  --merge A,B
```

### E) Prebuilt CG complex (`--complex-config`)

```bash
martinisurf \
  --complex-config input/complex_config.yaml \
  --surface-mode 4-1 \
  --lx 15 --ly 15 --dx 0.47 \
  --surface-bead P4 \
  --dist 10 \
  --substrate input/ETO.gro \
  --substrate-count 10 \
  --solvate --ionize --salt-conc 0.15
```

## 5) Flags you will use most often

| Flag | What it does | Typical value |
|---|---|---|
| `--pdb` | Input structure (local file, RCSB ID, or UniProt ID) | `1RJW` / `4C64.pdb` |
| `--moltype` | Molecule name for protein topology | `Protein` |
| `--dna` | Enables DNA workflow | no value |
| `--surface` | Reuses an existing surface file | `surface.gro` |
| `--surface-mode` | Builds 2-1 or 4-1 surface | `4-1` |
| `--lx --ly --dx` | Size and spacing for generated surface | `15 15 0.47` |
| `--anchor ...` | Residue-based anchor orientation | `--anchor 1 8 10 11` |
| `--linker` | Linker GRO file | `ALK.gro` |
| `--linker-group ...` | Residue groups for linker attachment | `--linker-group 1 1` |
| `--merge` | Chain merge during martinization | `A,B` / `A,B,C,D` |
| `--solvate` | Adds water with GROMACS | no value |
| `--ionize` | Adds ions (requires solvation) | no value |
| `--outdir` | Output folder | `Simulation_Files` |

## 6) Full flag reference (`martinisurf` main pipeline)

This section includes ALL flags from the main pipeline.

### Input and molecule

- `--pdb` (default: none): local `.pdb`, 4-character RCSB ID, or 6-character UniProt ID.
  - Required if you do NOT use `--complex-config`.
- `--complex-config` (default: none): YAML for pre-CG workflow (protein + cofactor), skipping martinization.
- `--moltype` (default: none): molecule name for protein topology.
- `--go` (default: false): enables Go model in `martinize2` (protein mode).
- `--ff` (default: `martini3001`): force field argument for `martinize2`.
- `--dna` (default: false): enables DNA mode (`martinize-dna.py`).
- `--dnatype` (default: `ds-stiff`): DNA type for `martinize-dna.py`.
- `--merge` (default: none, repeatable): chain merge (example: `A,B,C,D` or `all`).

### Martinization controls

- `--p {none,all,backbone}` (default: `backbone`): position restraint selection.
- `--pf` (default: `1000`): position restraint force constant.
- `--maxwarn` (default: `0`): max allowed warnings in `martinize2` before abort.
- `--dssp` (default: enabled): enables DSSP / secondary structure assignment in protein mode.
- `--no-dssp`: disables DSSP.
- `--elastic` (default: false): enables elastic network.
- `--ef` (default: `700`): elastic network force constant.
- `--go-eps` (default: none): Go epsilon.
- `--go-low` (default: none): Go minimum contact distance (nm).
- `--go-up` (default: none): Go maximum contact distance (nm).

### Surface

- `--surface` (default: none): existing surface `.gro` file.
- `--surface-mode {2-1,4-1}` (default: `2-1`): grid mode for generated surfaces.
- `--lx` (default: none): generated surface X size (nm).
- `--ly` (default: none): generated surface Y size (nm).
- `--dx` (default: `0.47`): bead spacing (2-1) or C-C parameter (4-1).
- `--surface-bead` (default: `C1`): bead type for generated surface.
- `--charge` (default: `0`): bead charge for generated surface topology.

### Classic orientation (anchors)

- `--anchor GROUP RESID [RESID ...]` (default: none, repeatable): defines anchor groups.
- `--dist` (default: `10.0` A): target anchor-to-surface distance.

### Linker orientation

- `--linker` (default: none): linker `.gro` file.
- `--linker-group GROUP RESID [RESID ...]` (default: none, repeatable): residue groups where linkers attach.
- `--linker-prot-dist` (default: auto): linker-to-biomolecule distance (A).
- `--linker-surf-dist` (default: auto): linker-to-surface distance (A).
- `--invert-linker` (default: false): reverses linker bead order.
- `--surface-linkers` (default: `0`): number of additional random surface linkers.

### Optional substrate

- `--substrate` (default: none): substrate `.gro` file for random insertion.
- `--substrate-itp` (default: inferred): substrate `.itp` file.
- `--substrate-count` (default: `0`): number of substrate molecules.

### Optional solvation and ionization

- `--solvate` (default: false): runs `gmx solvate`.
- `--ionize` (default: false): runs `gmx genion` after solvation.
- `--salt-conc` (default: `0.15` M): target salt concentration.
- `--water-gro` (default: none): custom solvent `.gro` for solvation.
- `--solvate-radius` (default: `0.21` nm): exclusion radius.
- `--solvate-surface-clearance` (default: `0.4` nm): removes water near surface plane.
- `--freeze-water-fraction` (default: `0.0`): fraction of `W` converted to `WF`.
  - Requires `--dna` and `--solvate`.
- `--freeze-water-seed` (default: `42`): seed parameter for freeze-water step.

### Output

- `--outdir` (default: `Simulation_Files`): output directory.

## 7) Advanced flags (internal modules)

These commands are usually called by MartiniSurf internally. Most new users do not need them directly.

### 7.1 `python -m martinisurf.surface_builder`

- `--mode {2-1,4-1}`
- `--bead`
- `--dx`
- `--lx`
- `--ly`
- `--lz`
- `--resname`
- `--output`
- `--charge`
- `--layers` (4-1)
- `--dist-z` (4-1)

### 7.2 `python -m martinisurf.system_tethered`

- `--surface`
- `--system`
- `--out`
- `--anchor`
- `--dist`
- `--linker-gro`
- `--linker-group`
- `--linker-prot-dist`
- `--linker-surf-dist`
- `--invert-linker`
- `--surface-linkers`
- `--surface-min-dist`
- `--dna-mode`

### 7.3 `python -m martinisurf.gromacs_inputs`

- `--moltype`
- `--outdir`
- `--anchor`
- `--linker-resid`
- `--use-linker`
- `--linker-resname`
- `--linker-size`
- `--linker-itp-name`
- `--linker-pull-init-prot`
- `--linker-pull-init-surf`
- `--go-model`
- `--cofactor-itp-name`
- `--cofactor-count`
- `--substrate-itp-name`
- `--substrate-count`

### 7.4 Legacy wrapper flags in `python -m martinisurf` (`orient` subcommand)

These exist for compatibility:
- `--surface`
- `--enzyme`
- `--out`
- `--dist`
- `--display`
- `--anchor`

Note: modern orientation module is `system_tethered` (uses `--system`, not `--enzyme`).

## 8) Output structure to check

Default output:

```text
Simulation_Files/
  0_topology/
    system.top
    system_res.top
    index.ndx
    system_itp/
  1_mdp/
  2_system/
    immobilized_system.gro
    system.gro
    (optional) final_system.gro
```

If you use `--solvate` and `--ionize`, check especially:
- `0_topology/system_final.top`
- `2_system/final_system.gro`
- `2_system/system_final.gro`

## 9) Common errors and direct fixes

- Error: `--pdb is required unless --complex-config is provided`
  - Fix: add `--pdb ...` or use `--complex-config ...`.

- Error: `When --surface is not provided, both --lx and --ly are required`
  - Fix: add `--lx` and `--ly`, or pass `--surface`.

- Error: `Linker mode requires at least one --linker-group`
  - Fix: add at least one `--linker-group`.

- Error: `--ionize requires --solvate`
  - Fix: enable `--solvate`.

- Error: `DNA mode requires python2.7`
  - Fix: install Python 2.7 or set `MARTINISURF_PYTHON2`.

- Error: missing linker `.itp`
  - Fix: place `linker.itp` next to `linker.gro` (same basename).

- Random substrate placement fails
  - Fix: reduce `--substrate-count` or use a larger box (`--lx`, `--ly`, `--dist`).

## 10) Practical checklist before long runs

- Visualize `2_system/immobilized_system.gro` before production simulation.
- Confirm `index.ndx` includes expected `Anchor_*` groups.
- Review `[ molecules ]` in `system.top` or `system_final.top`.
- If using linker/substrate/cofactor, verify their `.itp` files are in `0_topology/system_itp/`.
- Run a short test first (for example without `--ionize`) to validate geometry.

## 11) Short software evaluation

What is strong:
- End-to-end automation (structure to GROMACS-ready system).
- Supports protein, DNA, linker, and pre-CG complex workflows.
- Early validation for incompatible flag combinations.

What to keep in mind:
- Depends on external tools (`martinize2`, Python2, `gmx`).
- DNA path still depends on `martinize-dna.py` (Python2 ecosystem).
- `--freeze-water-seed` exists as a parameter, but current freeze-water selection is deterministic in the implementation.
