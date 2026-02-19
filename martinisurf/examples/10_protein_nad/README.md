# Pre-CG Protein + NAD + ETO (Solvate/Ionize) Example

Builds a surface-immobilized system from an already coarse-grained protein+cofactor complex, then:
- adds 10 random `ETO` molecules,
- solvates,
- ionizes.

Run:
```bash
bash run.sh
```

Input folder:
- `input/BsADH_NAD.gro`: pre-CG complex (protein + NAD).
- `input/BsADH_NAD.itp`: protein topology (`Active`).
- `input/NAD.itp`: cofactor topology (`NAD`).
- `input/go_*`: Gō includes.
- `input/ETO.gro` + `input/ETO.itp`: extra substrate molecule.
- `input/complex_config.yaml`: pre-CG workflow configuration.

Notes:
- `cofactor.count` is fixed to `4` in `complex_config.yaml`.
- Anchors are set as:
  - `--anchor 1 8 10 11`
  - `--anchor 2 1025 1027 1028`
- Uses generated surface (`--lx 20 --ly 20`), so no external `surface.gro` is required.
- Requires GROMACS (`gmx`) available in `PATH` for `--solvate --ionize`.
