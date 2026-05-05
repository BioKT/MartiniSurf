# MartiniSurf Quick Guide for New Users

This guide is for users with no prior Martini or coarse-grained simulation experience.

Goal:
- Know which command to run.
- Understand what each flag does.
- Avoid common setup errors.

## 1) 30-second rules

- Choose only ONE orientation mode:
  - `--anchor ...` (Not explicit Linker / anchor mode), or
  - `--linker ... --linker-group ...` (linker mode).
- `--ads-mode` is a specialized anchor workflow (not linker):
  uses anchor-based orientation but skips anchor pull/restraint topology and uses
  `minimization -> nvt -> npt -> production` (no `deposition` stage).
- In `--pdb` workflows, `--anchor` and `--linker-group` can use either global residue ids or `CHAIN RESID ...` syntax from the input PDB.
- In `--complex-config`, chain-based `protein.anchor_groups` also work if you provide `protein.reference_pdb`.
- For two-anchor systems, use `--balance-low-z` (and optional `--balance-low-z-fraction`) to flatten the lowest-Z protein face against the surface.
- In `--complex-config`, low-Z balancing is enabled by default (`protein.balance_low_z=true`), with default `protein.balance_low_z_fraction=0.2`.
- If you do NOT pass `--surface`, you must pass `--lx` and `--ly`.
- `--ionize` always requires `--solvate`.
- `--water-mix` works only with Martini 3 protein workflows and `--solvate`.
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
use protein mode + Not explicit Linker (`--anchor`).
- Protein with linker to surface:
use protein mode + `--linker` + `--linker-group`.
- Atomistic DNA PDB to surface:
use `--dna` + Not explicit Linker (`--anchor`) or explicit linker (`--linker`).
- Already CG protein+cofactor complex:
use `--complex-config`.

Recommended examples: `protein/04`, `dna/03`, and `protein/05` are the most complete solvated/ionized workflows:
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
- The pressure-coupled equilibration is done in `deposition`; there is no separate DNA `npt.mdp` stage in this example.

## 4) Copy-and-run commands

If you want the most complete and production-oriented setups, use examples `protein/04`, `dna/03`, and `protein/05` above.

### A) Example protein/04: Protein + Not explicit Linker (`--anchor`) + solvate + ionize

```bash
martinisurf \
  --pdb inputs/1RJW.pdb \
  --dssp \
  --go \
  --moltype Protein \
  --surface-mode 4-1 \
  --surface-bead P4 \
  --dx 0.53 \
  --lx 15 --ly 15 \
  --anchor A 8 10 11 \
  --anchor D 8 10 11 \
  --dist 1.0 \
  --solvate \
  --ionize \
  --salt-conc 0.15 \
  --merge A,B,C,D
```

### B) Example dna/03: DNA + linker + solvate + ionize + frozen water

```bash
martinisurf \
  --dna \
  --dnatype ds-stiff \
  --pdb inputs/4C64.pdb \
  --surface-mode 4-1 \
  --lx 10 \
  --ly 10 \
  --dx 0.27 \
  --surface-bead C1 \
  --linker inputs/ALK.gro \
  --linker-group A 1 \
  --solvate \
  --ionize \
  --salt-conc 0.15 \
  --freeze-water-fraction 0.10 \
  --freeze-water-seed 42 \
  --merge A,B
```

### C) Example protein/05: pre-CG protein+cofactor + substrate + solvate + ionize

```bash
martinisurf \
  --complex-config input/complex_config.yaml \
  --surface-mode 4-1 \
  --surface-bead P4 \
  --dx 0.53 \
  --lx 15 \
  --ly 15 \
  --dist 1.0 \
  --substrate input/ETO.gro \
  --substrate-count 10 \
  --solvate \
  --ionize \
  --salt-conc 0.15 \
  --outdir Simulation_Files
```

## 5) Flags you will use most often

| Flag | What it does | Typical value |
|---|---|---|
| `--pdb` | Input structure (local `.pdb`, `.cif`, `.mmcif`, RCSB ID, or UniProt ID) | `1RJW` / `4C64.pdb` / `model.cif` |
| `--moltype` | Molecule name for protein topology | `Protein` |
| `--dna` | Enables DNA workflow | no value |
| `--surface` | Reuses an existing surface file | `surface.gro` |
| `--surface-mode` | Builds `2-1`, `4-1`, `graphene`, `graphene-finite`, or `graphite` surfaces | `4-1` |
| `--lx --ly --dx` | Size and spacing for generated surface | `15 15 0.53` |
| `--surface-stacking` | Local multilayer stacking for `2-1`/`4-1` surfaces | `hcp` / `fcc` |
| `--anchor ...` | Not explicit Linker orientation (anchor mode) | `--anchor B 8 10 11` |
| `--ads-mode` | Anchor-like adsorption mode without anchor pull/restraint topology | no value |
| `--linker` | Linker GRO file | `ALK.gro` |
| `--linker-group ...` | Residue groups for linker attachment | `--linker-group A 1` |
| `--merge` | Chain merge during martinization | `A,B` / `A,B,C,D` |
| `--solvate` | Adds water with GROMACS | no value |
| `--ionize` | Adds ions (requires solvation) | no value |
| `--water-mix` | Martini 3 W/SW/TW water composition after solvation | `SW:0.10,TW:0.10` |
| `--outdir` | Output folder | `Simulation_Files` |

## 6) Full flag reference (`martinisurf` main pipeline)

This section includes ALL flags from the main pipeline.

### Input and molecule

- `--pdb` (default: none): local `.pdb`, `.cif`, or `.mmcif`, 4-character RCSB ID, or 6-character UniProt ID.
  - Local mmCIF/PDBx files are converted internally to PDB before the usual MartiniSurf cleaning and martinization steps.
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
- `--el`, `--eu`: elastic-network lower and upper distance cutoffs passed to `martinize2`.
- `--ermd`: minimum residue separation for elastic bonds.
- `--ea`, `--ep`: elastic-network decay factor and decay power.
- `--em`: minimum retained elastic-bond force constant.
- `--eb`: comma-separated bead names used for elastic bonds.
- `--eunit`: structural unit used to construct the elastic network.
- `--go-eps` (default: none): Go epsilon.
- `--go-low` (default: none): Go minimum contact distance (nm).
- `--go-up` (default: none): Go maximum contact distance (nm).
- `--go-res-dist`: minimum graph/residue distance for Go contacts.
- `--go-write-file [PATH]`: write the automatically calculated Go contact map.
- `--go-backbone`: backbone bead name used for Go virtual-site placement.
- `--go-atomname`: atom name used for Go virtual interaction sites.
- `--ss`: manual secondary-structure string passed to `martinize2`.
- `--collagen`: use `martinize2` collagen parameters.
- `--ed`: use extended-region dihedrals rather than elastic bonds.
- `--martinize-extra-args "..."`: advanced protein-mode passthrough to `martinize2`; MartiniSurf-managed input/output/topology flags are blocked.

### Surface

- `--surface` (default: none): existing surface `.gro` file.
- `--surface-mode {2-1,4-1,graphene,graphene-periodic,graphene-finite,graphite}` (default: `2-1`): mode for generated surfaces.
- `--surface-geometry {planar,3d}` (default: `planar`): orientation handling for the surface. Use `3d` for nanotubes and other non-planar surfaces.
- `--lx` (default: none): generated surface X size (nm).
- `--ly` (default: none): generated surface Y size (nm).
- `--dx` (default: `0.53`): in-plane nearest-neighbor lattice spacing in nm for local `2-1` / `4-1` surfaces. The generated lattice uses this value directly.
- `--surface-layers` (default: none): number of layers for local `2-1` / `4-1` modes.
- `--surface-stacking {hcp,fcc}` (default: `hcp`): local multilayer stacking for `2-1` / `4-1` surfaces. `hcp` uses ABAB stacking; `fcc` uses ABCABC stacking.
- `--surface-dist-z` (default: none): manual interlayer spacing in nm for local `2-1` / `4-1` surfaces. If omitted, MartiniSurf uses the close-packed HCP/FCC spacing `sqrt(2/3) * --dx`.
- `--surface-periodic-xy` / `--no-surface-periodic-xy` (default: enabled): connect local `2-1` / `4-1` surface lattice bonds across periodic X/Y boundaries.
- `--graphite-layers` (default: none): number of stacked graphene layers for `graphite`.
- `--graphite-spacing` (default: none): interlayer spacing in nm for `graphite`.
- `--surface-bead` (default: `C1`): bead type for generated surface.
- `--charge` (default: `0`): bead charge for generated surface topology.

Notes:
- `graphene` and `graphene-periodic` use the vendored Martini 3 periodic graphene generator.
- `graphene-finite` builds a finite graphene sheet.
- `graphite` builds a stacked graphitic slab and carries over the generated `posre_GRA.itp`.
- `graphene` and `graphite` are intended for Martini 3 / protein workflows.
- `planar` keeps the current surface orientation workflow unchanged.
- `3d` first lays the surface onto the `XY` plane using its principal axes and then places the biomolecule above the top external side.

Unit note:
- In the main `martinisurf` CLI, `--lx`, `--ly`, and `--dx` are provided in nm.

### Classic orientation (anchors)

- `--anchor GROUP_OR_CHAIN RESID [RESID ...]` (default: none, repeatable): defines anchor groups.
  - Legacy syntax: `--anchor 1 8 10 11`
  - Chain-based syntax for `--pdb` workflows: `--anchor D 8 10 11`
  - Chain-based groups are converted internally to global residue ids in appearance order (`Anchor_1`, `Anchor_2`, ...).
- `--dist` (default: `1.0` nm): target anchor-to-surface distance.
- `--ads-mode` (default: `false`): adsorption mode using anchor orientation without anchor pull/restraint topology.
  - MDP stages in this mode: `minimization`, `nvt`, `npt`, `production` (no `deposition`).
  - Incompatible with linker mode.
- `--balance-low-z` (default: `false`): in two-anchor mode, picks the roll angle that flattens the lowest-Z region.
- `--balance-low-z-fraction` (default: `0.2`): fraction (0,1] of lowest-Z beads used by `--balance-low-z`.

### Linker orientation

- `--linker` (default: none): linker `.gro` file.
- `--linker-group GROUP_OR_CHAIN RESID [RESID ...]` (default: none, repeatable): residue groups where linkers attach.
  - Legacy syntax: `--linker-group 1 8 10 11`
  - Chain-based syntax for `--pdb` workflows: `--linker-group B 8 10 11`
- `--linker-prot-dist` (default: auto): linker-to-biomolecule distance (nm).
  - Auto rule: `sigma` in DNA linker mode, `sigma * 1.2` in protein linker mode.
- `--linker-surf-dist` (default: auto): linker-to-surface distance (nm).
  - Auto rule: estimated from Martini bead-size sigma (`sigma * 1.2`) using linker-tail and surface bead classes.
- `--invert-linker` (default: false): reverses linker bead order.
- `--surface-linkers` (default: `0`): number of additional random surface linkers.
- DNA linker coupling behavior:
  - Uses bonded linker-DNA coupling in topology (no linker-DNA pull).
  - Bond target bead in DNA residue: `BB1`, else `BB2`, else `BB3`.
  - Adds linker-linker-DNA angle (`180`, force `20`).
  - Keeps linker-surface pull behavior unchanged.
  - Generated linker ITP includes linker-tail `position_restraints` under `#ifdef POSRES` for surface coupling control.

Notes:
- Chain-based syntax is resolved from the cleaned input PDB before martinization.
- In `--complex-config`, `protein.anchor_groups` can also use chain-based syntax if `protein.reference_pdb` points to the source PDB that matches the pre-CG complex residue order.
- In `--complex-config`, `protein.balance_low_z` defaults to `true`; optionally override with:
  - `protein.balance_low_z: false`
  - `protein.balance_low_z_fraction: 0.30`

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
- `--water-mix` (default: none): Martini 3 protein workflows only; converts final `W` waters into a requested `W`/`SW`/`TW` composition. Example: `--water-mix SW:0.10,TW:0.10` converts 10% to `SW`, 10% to `TW`, and keeps the remaining 80% as `W`.
- `--water-mix-seed` (default: `42`): seed parameter for deterministic W/SW/TW selection.
- `--freeze-water-fraction` (default: `0.0`): fraction of `W` converted to `WF`.
  - Requires `--dna` and `--solvate`.
- `--freeze-water-seed` (default: `42`): seed parameter for freeze-water step.

### Output

- `--outdir` (default: `Simulation_Files`): output directory.

## 7) Advanced flags (internal modules)

These commands are usually called by MartiniSurf internally. Most new users do not need them directly.

### 7.1 `python -m martinisurf.surface_builder`

- `--mode {2-1,4-1,graphene,graphene-periodic,graphene-finite,graphite,cnt-m2,cnt-m3}`
- `--bead`
- `--dx`
- `--lx`
- `--ly`
- `--lz`
- `--resname`
- `--output`
- `--charge`
- `--layers` for local multilayer `2-1` / `4-1` surfaces
- `--stacking {hcp,fcc}` for local multilayer `2-1` / `4-1` surfaces
- `--dist-z` manual interlayer spacing for local multilayer `2-1` / `4-1` surfaces
- `--periodic-xy` to connect local lattice bonds across periodic X/Y boundaries
- `--graphite-layers`
- `--graphite-spacing`
- `--cnt-numrings`
- `--cnt-ringsize`
- `--cnt-bondlength`
- `--cnt-bondforce`
- `--cnt-angleforce`
- `--cnt-beadtype`
- `--cnt-functype`
- `--cnt-func-begin`
- `--cnt-func-end`
- `--cnt-base36`

Notes:
- `cnt-m2` wraps the vendored `cnt-martini` default Martini 2 preset.
- `cnt-m3` uses the same CNT generator with Martini 3-style defaults: `0.41 nm` bond length, `9` beads per ring, and `SC5` carbon beads.

### 7.2 `python -m martinisurf.system_tethered`

- `--surface`
- `--system`
- `--out`
- `--anchor`
- `--anchor-landmark-mode`
- `--dist`
- `--reference-exclude-resname`
- `--linker-gro`
- `--linker-group`
- `--linker-prot-dist`
- `--linker-surf-dist`
- `--invert-linker`
- `--surface-linkers`
- `--surface-min-dist`
- `--dna-mode`
- `--min-reference-z-dist`
- `--balance-low-z`
- `--balance-low-z-fraction`

### 7.3 `python -m martinisurf.gromacs_inputs`

- `--moltype`
- `--outdir`
- `--anchor`
- `--linker-resid`
- `--use-linker`
- `--linker-resname`
- `--linker-size`
- `--linker-itp-name`
- `--linker-pull-init-prot` (optional): explicit initial distance (nm) for linker-head to biomolecule pull coordinates.
- `--linker-pull-init-surf` (optional): explicit initial distance (nm) for linker-tail to surface pull coordinates.
- `--go-model`
- `--ads-mode`
- `--cofactor-itp-name`
- `--cofactor-count`
- `--substrate-itp-name`
- `--substrate-moltype`
- `--substrate-count`

MDP pull note:
- In generated pull coordinates, MartiniSurf writes either `pull-coordX_start = yes` or `pull-coordX-init = ...`, never both for the same coordinate.
- `--linker-pull-init-prot` and `--linker-pull-init-surf` in `gromacs_inputs` are in nm (internal module API).
- In DNA linker mode, linker-DNA is handled with bonded terms in topology instead of pull coordinates.

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
  provenance.json
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

`provenance.json` records the command-line arguments, workflow mode, key input paths,
selected model parameters, dependency versions/paths when detectable, and the boundary
between MartiniSurf-managed steps and user/external parametrizations.

If you use `--solvate` and `--ionize`, check especially:
- `0_topology/system_final.top`
- `0_topology/system_final_res.top`
- `2_system/final_system.gro`
- `2_system/system_final.gro`

Workflow note for example `work_flow_gromacs.sh` scripts:
- `minimization`, `nvt`, `npt`, and `deposition` use `system_final.top` (or base non-restrained topology).
- `production` uses `system_final_res.top` when present (anchor/linker focused restraints).

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
- `--freeze-water-seed` is active in the current implementation and controls the stochastic water-to-`WF` selection when `--freeze-water-fraction` is used.

## 12) Colab behavior notes

- In `MartiniSurf_Protein.ipynb` and `MartiniSurf_DNA.ipynb`, Step `6B` runs `grompp` with `-maxwarn 3`.
- In `MartiniSurf_DNA.ipynb`, Step `6B` uses the DNA short-MD protocol `nvt -> deposition (NPT) -> production (NVT)`.
- Step `6C - View MD Result` renders the stages that are actually available for the selected notebook workflow.
