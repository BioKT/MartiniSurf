# Protein Linker + Surface Decoration Example

Builds a protein-on-surface system using linker mode and decorates the surface with random extra linkers.

Run:
```bash
bash run.sh
```

Notes:
- Uses `--go --dssp` explicitly for protein martinization.
- Uses `inputs/EPOXY.gro` as linker.
- `inputs/EPOXY.itp` must match linker basename (`EPOXY.gro` -> `EPOXY.itp`).
- Uses two linker groups as attachment targets.
- `--surface-linkers 12` adds 12 random linkers above the surface.
