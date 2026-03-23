# Protein Anchor + Solvate + Ionize Example

Builds a protein-on-surface system in anchor mode, then runs optional GROMACS solvation and ionization.

Run:
```bash
bash run.sh
```

Notes:
- Uses `--go --dssp` explicitly for protein martinization.
- Uses chain-local anchor syntax in `run.sh` (`--anchor A 8 10 11` and `--anchor D 8 10 11`).
- Requires GROMACS (`gmx`) installed and available in `PATH`.
- Generates `0_topology/system_final.top` and `2_system/final_system.gro`.
