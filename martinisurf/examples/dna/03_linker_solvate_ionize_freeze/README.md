# DNA Linker + Solvate + Ionize + 10% Frozen Water Example

Builds a DNA-on-surface system in linker mode, runs optional GROMACS solvation/ionization, then uses MartiniSurf DNA water freezing to convert 10% of water molecules from `W` to `WF`.

Run:
```bash
bash run.sh
```

What this example configures in `martinisurf`:
- Generated surface: `--surface-mode 4-1 --lx 10 --ly 10 --dx 0.27 --surface-bead C1`
- Linker from inputs: `--linker inputs/ALK.gro --linker-group A 1` (uses `inputs/ALK.itp`)
- `--freeze-water-fraction 0.10`
- `--freeze-water-seed 42`

DNA linker coupling details:
- Linker-DNA uses bonded coupling in topology (`bond + angle`) instead of pull.
- DNA target bead priority per selected residue: `BB1`, else `BB2`, else `BB3`.
- Linker-surface pull remains active.

Generated outputs are updated automatically:
- `Simulation_Files/2_system/final_system.gro`
- `Simulation_Files/2_system/system_final.gro`
- `Simulation_Files/0_topology/system_final.top`
- `Simulation_Files/0_topology/system_final_res.top`

Workflow note:
- `work_flow_gromacs.sh` uses `system_final.top` for equilibration stages and `system_final_res.top` for production when present.

Notes:
- Assumes Martini solvent files already define `W` and `WF` behavior in your FF stack.
- Requires GROMACS (`gmx`) installed and available in `PATH`.
