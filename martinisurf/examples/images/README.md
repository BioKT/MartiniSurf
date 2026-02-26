# Example Renders (VMD + Tachyon)

Generate one final render per numbered example folder (`01_*`, `02_*`, ...):

```bash
bash martinisurf/examples/images/run_render.sh
```

Output files are written to:
- `martinisurf/examples/images/*_top.png`
- `martinisurf/examples/images/*_side.png`
- `martinisurf/examples/images/*_final.png` (alias to side view for compatibility)

Render rules:
- Projection: orthographic
- View: ZX-style
- Ray tracing: `TachyonInternal` with AO/shadows enabled
- Protein: only `name BB`
- Protein anchors (example 05 residues `8 10 11 1025 1027 1028`): highlighted in red
- Solvated systems: water hidden (`resname W WF SOL`)
