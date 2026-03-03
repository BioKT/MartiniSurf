# DNA Linker Example

Builds a DNA-on-surface system using linker mode.

Run:
```bash
bash run.sh
```

Notes:
- `inputs/linker.itp` must be present with the same basename as `linker.gro`.
- `run.sh` uses chain-local residue syntax (`--linker-group A 1`).
- Add more `--linker-group` entries in `run.sh` for multiple linkers.
