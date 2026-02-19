#!/usr/bin/env python3
"""
Surface Builder Module
Generates 2D surfaces for Martini simulations using two mapping modes:
- 2-1: Standard hexagonal mapping.
- 4-1: Honeycomb Carbon mapping (atomistic-like lattice).
"""

import argparse
import math
import os
import shutil
from pathlib import Path
from typing import Iterable, List, Tuple

def main(argv: Iterable[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Martini Surface Generator (2-1 and 4-1 mappings)")

    # Base parameters
    parser.add_argument("--mode", type=str, choices=["2-1", "4-1"], default="2-1", 
                        help="Mapping mode: 2-1 (Standard hex) or 4-1 (Honeycomb Carbon)")
    parser.add_argument("--bead", type=str, default="P4", help="Martini bead type")
    parser.add_argument("--dx", type=float, default=0.47, help="Bead spacing (2-1) or C-C distance (4-1)")
    parser.add_argument("--lx", type=float, required=True, help="Surface length in X (nm)")
    parser.add_argument("--ly", type=float, required=True, help="Surface length in Y (nm)")
    parser.add_argument("--lz", type=float, default=10.0, help="Box height in Z (nm)")
    parser.add_argument("--resname", type=str, default="SRF", help="Residue name")
    parser.add_argument("--output", type=str, default="surface", help="Output basename")
    parser.add_argument("--charge", type=float, default=0.0, help="Charge per bead")

    # Honeycomb Carbon (4-1) specific parameters
    parser.add_argument("--layers", type=int, default=1, help="Number of layers (4-1 mode only)")
    parser.add_argument("--dist-z", type=float, default=0.335, help="Interlayer spacing in nm")

    args = parser.parse_args(list(argv) if argv is not None else None)

    outdir = os.path.dirname(args.output) or "."
    basename = os.path.basename(args.output)
    os.makedirs(outdir, exist_ok=True)

    atoms: List[Tuple[float, float, float]] = []
    final_lx, final_ly = args.lx, args.ly

    # =========================================================
    # MODE 2-1: STANDARD HEXAGONAL MAPPING
    # =========================================================
    if args.mode == "2-1":
        scale = args.dx / 0.142
        a = 0.246 * scale
        atoms_unit = [
            (0.0, 0.0, 0.0), (a, 0.0, 0.0), (2 * a, 0.0, 0.0),
            (0.5 * a, (math.sqrt(3) / 2) * a, 0.0),
            (1.5 * a, (math.sqrt(3) / 2) * a, 0.0),
            (2.5 * a, (math.sqrt(3) / 2) * a, 0.0),
        ]
        lx_cell, ly_cell = 3.0 * a, math.sqrt(3) * a
        nx, ny = max(1, round(args.lx / lx_cell)), max(1, round(args.ly / ly_cell))
        final_lx, final_ly = nx * lx_cell, ny * ly_cell

        for i in range(nx):
            for j in range(ny):
                for x, y, z in atoms_unit:
                    atoms.append((x + i * lx_cell, y + j * ly_cell, 3.0))

    # =========================================================
    # MODE 4-1: HONEYCOMB CARBON MAPPING
    # =========================================================
    else:
        d_cc = args.dx # In 4-1 mode, dx is treated as C-C distance
        lx_cell, ly_cell = math.sqrt(3) * d_cc, 3 * d_cc
        nx, ny = max(1, round(args.lx / lx_cell)), max(1, round(args.ly / ly_cell))
        final_lx, final_ly = nx * lx_cell, ny * ly_cell

        # Unit cell atoms for Honeycomb lattice
        base_unit = [
            [0.0, 0.0, 0.0], [lx_cell/2, d_cc/2, 0.0],
            [lx_cell/2, 1.5*d_cc, 0.0], [0.0, 2*d_cc, 0.0]
        ]
        
        for layer in range(args.layers):
            # AB (Bernal) Stacking shift
            shift_x = (lx_cell / 2) if (layer % 2 != 0) else 0.0
            shift_y = d_cc if (layer % 2 != 0) else 0.0
            z_pos = 3.0 + (layer * args.dist_z)
            
            for i in range(nx):
                for j in range(ny):
                    for ax, ay, az in base_unit:
                        atoms.append((ax + i * lx_cell + shift_x, ay + j * ly_cell + shift_y, az + z_pos))

    # ---------------------------------------------------------
    # File Writing
    # ---------------------------------------------------------
    gro_path = Path(outdir, f"{basename}.gro")
    with gro_path.open("w") as fgro:
        fgro.write(f"Surface {args.mode} | {args.lx}x{args.ly} nm | Layers: {args.layers}\n")
        fgro.write(f"{len(atoms):5d}\n")
        for i, (x, y, z) in enumerate(atoms, start=1):
            fgro.write(f"{1:5d}{args.resname:<4}{args.bead:>4}{i:6d}{x:8.3f}{y:8.3f}{z:8.3f}\n")
        fgro.write(f"{final_lx:10.5f}{final_ly:10.5f}{args.lz:10.5f}\n")

    itp_path = Path(outdir, f"{basename}.itp")
    with itp_path.open("w") as fitp:
        fitp.write(";;;;;; Minimal surface topology\n\n[ moleculetype ]\n; molname nrexcl\n")
        fitp.write(f"  {args.resname}        1\n\n[ atoms ]\n")
        fitp.write(f"  1   {args.bead:<6}   1   {args.resname:<4}   C     1     {args.charge:.3f}\n")

    # Backward-compatible copy for standalone tooling/tests.
    active_dst = resolve_activeitp_destination(outdir)
    active_dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(itp_path, active_dst)

    print(f"✔ {args.mode} Surface ({'Honeycomb Carbon' if args.mode == '4-1' else 'Hex Mapping'}) generated with {len(atoms)} beads.")

def resolve_activeitp_destination(outdir: str) -> Path:
    outdir_path = Path(outdir).resolve()
    default_dst = outdir_path / "system_itp" / "surface.itp"

    # Prefer centralized topology includes when called from MartiniSurf pipelines.
    for base in [outdir_path, *outdir_path.parents]:
        topology_dir = base / "0_topology" / "system_itp"
        if topology_dir.exists():
            return topology_dir / "surface.itp"

    return default_dst

if __name__ == "__main__":
    main()
