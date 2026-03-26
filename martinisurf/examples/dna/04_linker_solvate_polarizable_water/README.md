# DNA Linker + Solvate + Polarizable Water Example

Builds a DNA-on-surface system in linker mode, runs GROMACS solvation + ionization, and uses Martini 2 polarizable water (`PW`) instead of standard `W` water.

Run:
```bash
bash run.sh
```

What this example configures in `martinisurf`:
- Generated surface: `--surface-mode 4-1 --lx 10 --ly 10 --dx 0.27 --surface-bead C1`
- Linker from inputs: `--linker inputs/ALK.gro --linker-group A 1` (uses `inputs/ALK.itp`)
- Polarizable water mode: `--polarizable-water`
- Solvation: `--solvate`
- Ionization: `--ionize --salt-conc 0.15`
- Surface clearance: uses the default DNA value (`0.4`)

Polarizable-water specifics:
- Solvent template: `polarize-water.gro`
- DNA force-field include: `#include "martini_v2.1P-dna.itp"`
- Solvent molecule in topology: `PW`
- DNA `.mdp` files are rewritten to use the polarizable-water electrostatics/VdW and bond settings.

DNA linker coupling details:
- Linker-DNA uses bonded coupling in topology (`bond + angle`) instead of pull.
- DNA target bead priority per selected residue: `BB1`, else `BB2`, else `BB3`.
- Linker-surface pull remains active.

Generated outputs are updated automatically:
- `Simulation_Files/2_system/solvated_system.gro`
- `Simulation_Files/2_system/final_system.gro`
- `Simulation_Files/2_system/system_final.gro`
- `Simulation_Files/0_topology/system_final.top`
- `Simulation_Files/0_topology/system_final_res.top`

Compatibility note:
- The generated polarizable-water `.mdp` files use a modern GROMACS/Verlet-compatible approximation of the legacy Martini 2 PW setup.
- The example now includes ionization and produces `PW`, `NA`, and `CL` in the final topology/system.

Notes:
- This example is intended for Martini 2 DNA with polarizable water only.
- `--freeze-water-fraction` is intentionally not used here because it is incompatible with `--polarizable-water`.
- Requires GROMACS (`gmx`) installed and available in `PATH`.
