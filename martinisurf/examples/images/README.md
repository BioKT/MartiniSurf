# Example Renders (VMD + Tachyon)

Generate one final render per numbered example folder (`01_*`, `02_*`, ...):

```bash
bash martinisurf/examples/images/run_render.sh
```

Output files are written to:
- `martinisurf/examples/images/*_final.png` (or `.tga` if PNG conversion is unavailable)

Render rules:
- Projection: orthographic
- View: ZX-style
- Ray tracing: `TachyonInternal` with AO/shadows enabled
- Protein: only `name BB`
- Solvated systems: water hidden (`resname W WF SOL`)
