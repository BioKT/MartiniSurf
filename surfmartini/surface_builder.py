#!/usr/bin/env python3
"""
Generate a hexagonal surface and a minimal GROMACS topology (.itp).

This module creates a 2D hexagonal lattice of beads to represent
a flat surface for Martini or GōMartini simulations.

It outputs:
- A .gro coordinate file
- A minimal .itp topology file
- (Standalone mode only) Automatic copy of surface.itp into ActiveITP/
"""

import argparse
import math
import os
import shutil
from pathlib import Path
from typing import Iterable, List, Tuple


# ======================================================================
# Entry point
# ======================================================================
def main(argv: Iterable[str] | None = None) -> None:
    """CLI entry point for hexagonal surface generation."""
    parser = argparse.ArgumentParser(
        description="Generate a hexagonal surface and minimal surface.itp"
    )

    parser.add_argument("--bead", type=str, default="P4")
    parser.add_argument("--dx", type=float, default=0.47)
    parser.add_argument("--lx", type=float, required=True)
    parser.add_argument("--ly", type=float, required=True)
    parser.add_argument("--lz", type=float, default=5.0)
    parser.add_argument("--resname", type=str, default="SRF")

    parser.add_argument(
        "--output",
        type=str,
        default="surface",
        help="Basename or full path for the output files",
    )

    parser.add_argument("--charge", type=float, default=0.0)

    args = parser.parse_args(list(argv) if argv is not None else None)

    # ---------------------------------------------------------
    # Resolve output directory and name
    # ---------------------------------------------------------
    outdir = os.path.dirname(args.output) or "."
    basename = os.path.basename(args.output)
    os.makedirs(outdir, exist_ok=True)

    # ---------------------------------------------------------
    # Hexagonal lattice parameters
    # ---------------------------------------------------------
    scale = args.dx / 0.142
    a = 0.246 * scale  # primitive constant

    atoms_unit: List[Tuple[float, float, float]] = [
        (0.0, 0.0, 0.0),
        (a, 0.0, 0.0),
        (2 * a, 0.0, 0.0),
        (0.5 * a, (math.sqrt(3) / 2) * a, 0.0),
        (1.5 * a, (math.sqrt(3) / 2) * a, 0.0),
        (2.5 * a, (math.sqrt(3) / 2) * a, 0.0),
    ]

    lx_cell = 3.0 * a
    ly_cell = math.sqrt(3) * a

    nx = max(1, round(args.lx / lx_cell))
    ny = max(1, round(args.ly / ly_cell))

    lx = nx * lx_cell
    ly = ny * ly_cell

    print(f"• Building hex surface: {lx:.3f} × {ly:.3f} × {args.lz:.3f} nm")

    # ---------------------------------------------------------
    # Build atom coordinates
    # ---------------------------------------------------------
    atoms: List[Tuple[float, float, float]] = []
    for i in range(nx):
        for j in range(ny):
            dx = i * lx_cell
            dy = j * ly_cell
            for x, y, z in atoms_unit:
                atoms.append((x + dx, y + dy, z))

    # ---------------------------------------------------------
    # Write GRO file
    # ---------------------------------------------------------
    gro_path = Path(outdir, f"{basename}.gro")

    with gro_path.open("w") as fgro:
        fgro.write(f"Hex surface {args.lx}x{args.ly} nm (d={args.dx:.2f} nm)\n")
        fgro.write(f"{len(atoms):5d}\n")

        for i, (x, y, z) in enumerate(atoms, start=1):
            fgro.write(
                f"{i:5d}{args.resname:<4}{args.bead:>4}{i:6d}"
                f"{x:8.3f}{y:8.3f}{z:8.3f}\n"
            )

        fgro.write(f"{lx:10.5f}{ly:10.5f}{args.lz:10.5f}\n")

    print(f"✔ Wrote {gro_path} ({len(atoms)} beads)")

    # ---------------------------------------------------------
    # Write ITP file
    # ---------------------------------------------------------
    itp_path = Path(outdir, f"{basename}.itp")

    with itp_path.open("w") as fitp:
        fitp.write(";;;;;; Minimal surface topology\n\n")
        fitp.write("[ moleculetype ]\n; molname   nrexcl\n")
        fitp.write(f"  {args.resname}        1\n\n")
        fitp.write("[ atoms ]\n")
        fitp.write("; id  type  resnr  residu  atom  cgnr  charge\n")
        fitp.write(
            f"  1   {args.bead:<6}   1   {args.resname:<4}   C     1     {args.charge:.3f}\n"
        )

    print(f"✔ Wrote {itp_path}")

    # ---------------------------------------------------------
    # Copy surface.itp (ONLY in standalone mode)
    # ---------------------------------------------------------
    final_dst = resolve_activeitp_destination(outdir)

    # If path indicates pipeline mode → skip
    if "2_system" in str(final_dst):
        print("• Pipeline mode detected — skipping auto-copy of surface.itp")
    else:
        final_dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(itp_path, final_dst)
        print(f"✔ Copied surface.itp → {final_dst}")


# ======================================================================
# Helper: detect correct destination for surface.itp
# ======================================================================
def resolve_activeitp_destination(outdir: str) -> Path:
    """
    Determine correct surface.itp destination:
    - Standalone mode: outdir/ActiveITP/surface.itp
    - Pipeline mode:   SKIP (return a path inside 2_system but caller won't use it)
    """
    outdir_path = Path(outdir)
    default_dst = outdir_path / "ActiveITP" / "surface.itp"

    # Pipeline detection
    if "Simulation" in str(outdir_path):
        parts = list(outdir_path.parts)
        if "Simulation" in parts:
            sim_root = Path(*parts[: parts.index("Simulation") + 1])
            topology_dir = sim_root / "0_topology" / "ActiveITP"
            if topology_dir.exists():
                return topology_dir / "surface.itp"

    return default_dst


# ======================================================================
if __name__ == "__main__":
    main()
