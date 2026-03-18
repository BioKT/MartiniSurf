# Protein Immobilization on a CNT Example

Builds a protein-on-nanotube system in a single `martinisurf` run. The main pipeline generates the CNT automatically from `--surface-mode cnt` and then switches to `3d` surface orientation for the nanotube.

Run:
```bash
bash run.sh
```

Workflow:
1. Use `--surface-mode cnt` in the main pipeline.
2. MartiniSurf resolves that to the CNT preset that matches the workflow (`cnt-m3` for protein, `cnt-m2` for DNA).
3. The CNT is generated automatically and oriented with the `3d` surface workflow.

Key options:
- `--surface-mode cnt`: auto-selects the CNT preset in the main pipeline.
- `--cnt-numrings 24`: increases the nanotube length.
- `--cnt-ringsize`: changes the nanotube diameter.
- `--dist 1.0`: target anchor-to-surface distance (nm).

Notes:
- Protein workflows default to `cnt-m3` (`0.41 nm`, `9` beads per ring, `SC5` beads).
- DNA workflows default to `cnt-m2` (`0.47 nm`, `8` beads per ring, `CNP` beads).
- The CNT builder still requires `--lx` and `--ly` at the CLI level, but the nanotube geometry is controlled mainly by `--cnt-*`.
