# Overview

MartiniSurf builds complete GROMACS-ready systems for:

- Protein-surface workflows (Martini 3)
- DNA-surface workflows (Martini2 DNA)
- Linker and non-linker orientation workflows
- Pre-CG complex workflows (`--complex-config`)

## Main capabilities

- Coarse graining with `martinize2` (protein) and `martinize-dna.py` (DNA)
- Generated or user-provided surfaces
- Orientation by residue anchors (Not explicit Linker mode)
- Explicit linker mode
- Automatic topology assembly
- Optional solvation, ionization, and DNA water freezing

## Output structure

Typical output folder:

```text
Simulation_Files/
  0_topology/
  1_mdp/
  2_system/
```
