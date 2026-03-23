# Protein Adsorption Example (`--ads-mode`)

Builds a protein-on-surface system using adsorption mode (`--ads-mode`).

Run:
```bash
bash run.sh
```

Key options:
- `--dssp`: enables DSSP-assisted secondary structure assignment
- `--anchor`: two anchor groups selected by chain-local residue ids (`A 8 10 11` and `D 8 10 11`)
- `--ads-mode`: adsorption workflow without anchor pull/restraint topology
- `--dist 1.0`: target reference-to-surface distance (nm)
