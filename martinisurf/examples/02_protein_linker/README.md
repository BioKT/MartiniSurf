# Protein Linker Example

Builds a protein-on-surface system using linker mode.

Run:
```bash
bash run.sh
```

Notes:
- Uses `--go --dssp` explicitly for protein martinization.
- `inputs/linker.itp` must match linker basename (`linker.gro` -> `linker.itp`).
- Uses two linker groups as attachment targets.
