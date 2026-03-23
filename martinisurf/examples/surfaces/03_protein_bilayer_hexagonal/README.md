# Protein Immobilization on a Two-Layer Hexagonal Surface

Builds the same protein used in `../01_protein_cnt_3d` and `../02_protein_graphene_resizable` (`2A3D`), but on a two-layer local hexagonal surface.

Run:
```bash
bash run.sh
```

Key options:
- `--surface-mode 2-1`: uses the local hexagonal lattice
- `--surface-layers 2`: builds two stacked layers
- `--surface-dist-z 0.382`: spacing between layers
- `--surface-bead P4 C1`: first layer uses `P4`, second layer uses `C1`

Notes:
- You can change `--surface-bead P4 C1` to any sequence of bead types. They are assigned by layer in cyclic order.
- You can change `--lx` and `--ly` in `run.sh` to resize the surface.
- If `2A3D` is not cached locally, the run requires internet access to download it from RCSB.
