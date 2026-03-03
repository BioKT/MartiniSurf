# Protein Anchor Example

Builds a protein-on-surface system using classical anchor mode.

Run:
```bash
bash run.sh
```

Key options:
- `--go`: enables Go model (protein only)
- `--dssp`: enables DSSP-assisted secondary structure assignment
- `--anchor`: two anchor groups selected by chain-local residue ids (`A 8 10 11` and `D 8 10 11`)
- `--dist 10`: target anchor-to-surface distance
