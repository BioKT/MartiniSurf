# Protein Anchor + Solvate + Ionize Example

Builds a protein-on-surface system in anchor mode, then runs optional GROMACS solvation and ionization.

Run:
```bash
bash run.sh
```

Notes:
- Uses `--go --dssp` explicitly for protein martinization.
- Requires GROMACS (`gmx`) installed and available in `PATH`.
- Generates `0_topology/system_final.top` and `2_system/final_system.gro`.
