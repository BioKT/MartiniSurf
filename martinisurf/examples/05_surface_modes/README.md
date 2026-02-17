# Surface Modes Example (2-1 and 4-1)

Generates two standalone surfaces so you can inspect and compare both mapping modes:

- `2-1`: standard hexagonal mapping
- `4-1`: honeycomb carbon mapping

Run:
```bash
bash run.sh
```

Outputs are written into `outputs/`:
- `outputs/surface_2to1.gro`
- `outputs/surface_2to1.itp`
- `outputs/surface_4to1.gro`
- `outputs/surface_4to1.itp`

Quick checks:
- Open both `.gro` files and compare atom count (line 2)
- Compare final box line (last line)
- In `4-1` mode, increase `--layers` to generate multilayer surfaces
