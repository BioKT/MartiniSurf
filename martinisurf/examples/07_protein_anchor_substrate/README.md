# Protein Anchor + Substrate Example

Builds a protein-on-surface system in anchor mode and adds random coarse-grained substrate molecules in the box.

Run:
```bash
bash run.sh
```

Notes:
- Uses `--go --dssp` explicitly for protein martinization.
- `--substrate` infers `substrate.itp` from the same basename.
- `--substrate-count 20` inserts 20 random substrate molecules.
