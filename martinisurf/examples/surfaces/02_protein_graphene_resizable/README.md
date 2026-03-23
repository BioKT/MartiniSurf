# Protein Immobilization on Resizable Graphene

Builds the same protein used in `../01_protein_cnt_3d`, but immobilized on a graphene sheet instead of a nanotube.

Run:
```bash
bash run.sh
```

Resize the graphene by editing these lines in `run.sh`:
```bash
  --lx 20 \
  --ly 12 \
```

Key options:
- `--lx`: graphene length in X (nm)
- `--ly`: graphene length in Y (nm)
- `--dist`: anchor-to-surface distance (nm)
- `--outdir`: output folder

Workflow:
1. Downloads or reads the same protein used in surface example 01 (`2A3D`).
2. Generates a graphene surface directly inside the main `martinisurf` pipeline.
3. Immobilizes the protein on the planar graphene sheet with anchor mode.

Notes:
- Default graphene size is `15 x 15 nm`.
- This example keeps the same anchor selection as surface example 01: `--anchor A 1`.
- If you want a larger or smaller adsorption area, change `--lx` and `--ly` directly in `run.sh`.
- If `2A3D` is not cached locally, the run requires internet access to download it from RCSB.
