#!/usr/bin/env python3
"""
Surface Builder Module
Generates 2D surfaces for Martini simulations using generic lattice modes
and vendored Martini 3 carbon surface builders:
- 2-1: Standard hexagonal mapping.
- 4-1: Honeycomb Carbon mapping (atomistic-like lattice).
- graphene / graphene-finite: Martini 3 graphene sheets.
- graphite: Stacked graphitic slab.
"""

import argparse
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable, List, Tuple

LOCAL_SURFACE_MODES = {"2-1", "4-1"}
GRAPHENE_SURFACE_MODES = {
    "graphene": "graphene",
    "graphene-periodic": "graphene",
    "graphene-finite": "graphene-finite",
    "graphite": "graphite",
}
CNT_SURFACE_MODES = {
    "cnt-m2": "cnt-m2",
    "cnt-martini2": "cnt-m2",
    "cnt-m3": "cnt-m3",
    "cnt-martini3": "cnt-m3",
}
VALID_SURFACE_MODES = (
    "2-1",
    "4-1",
    "graphene",
    "graphene-periodic",
    "graphene-finite",
    "graphite",
    "cnt-m2",
    "cnt-martini2",
    "cnt-m3",
    "cnt-martini3",
)


def normalize_surface_mode(mode: str) -> str:
    clean_mode = (mode or "").strip().lower()
    if clean_mode in LOCAL_SURFACE_MODES:
        return clean_mode
    if clean_mode in GRAPHENE_SURFACE_MODES:
        return GRAPHENE_SURFACE_MODES[clean_mode]
    return CNT_SURFACE_MODES.get(clean_mode, "")


def repo_graphene_root() -> Path:
    return Path(__file__).resolve().parent / "surfaces_generator" / "Martini3-Graphene"


def repo_cnt_root() -> Path:
    return Path(__file__).resolve().parent / "surfaces_generator" / "cnt-martini"


def copy_support_files(extra_files: list[Path], destination_dir: Path) -> None:
    destination_dir.mkdir(parents=True, exist_ok=True)
    for extra in extra_files:
        if extra.exists():
            shutil.copy(extra, destination_dir / extra.name)


def rewrite_cnt_posres_include(itp_path: Path, posres_name: str) -> None:
    lines = []
    replaced = False
    for raw in itp_path.read_text().splitlines():
        if raw.strip().startswith('#include "') and raw.strip().endswith('-posres.itp"'):
            lines.append(f'#include "{posres_name}"')
            replaced = True
        else:
            lines.append(raw)
    if replaced:
        itp_path.write_text("\n".join(lines).rstrip() + "\n")


def normalize_layer_beads(bead_values: str | list[str] | tuple[str, ...] | None) -> list[str]:
    if bead_values is None:
        return ["P4"]
    if isinstance(bead_values, str):
        raw_values = [bead_values]
    else:
        raw_values = list(bead_values)
    cleaned = [str(value).strip() for value in raw_values if str(value).strip()]
    return cleaned or ["P4"]


def bead_for_layer(layer_beads: list[str], layer_index: int) -> str:
    return layer_beads[layer_index % len(layer_beads)]


def graphite_ab_shift(layer_index: int, x_shift: float, y_shift: float) -> tuple[float, float]:
    if layer_index % 2 == 0:
        return 0.0, 0.0
    return x_shift, y_shift


def write_local_surface_itp(
    itp_path: Path,
    resname: str,
    atoms: list[tuple[float, float, float, str]],
    charge: float,
) -> None:
    itp_lines = [
        ";;;;;; Minimal surface topology",
        "",
        "[ moleculetype ]",
        "; molname nrexcl",
        f"  {resname}        1",
        "",
        "[ atoms ]",
    ]

    unique_beads = {bead for _, _, _, bead in atoms}
    if len(unique_beads) == 1:
        bead = atoms[0][3] if atoms else "C1"
        itp_lines.append(
            f"  1   {bead:<6}   1   {resname:<4}   {bead:<5} 1     {charge:.3f}"
        )
    else:
        itp_lines.append("; nr type resnr residue atom cgnr charge")
        for idx, (_, _, _, bead) in enumerate(atoms, start=1):
            atom_name = f"B{idx % 100000}"
            itp_lines.append(
                f"  {idx:<5d}{bead:<7}1   {resname:<4}   {atom_name:<5} {idx:<5d} {charge:.3f}"
            )

    itp_path.write_text("\n".join(itp_lines).rstrip() + "\n")


def run_cnt_generator(
    mode: str,
    outdir: Path,
    basename: str,
    cnt_numrings: int | None,
    cnt_ringsize: int | None,
    cnt_bondlength: float | None,
    cnt_bondforce: float | None,
    cnt_angleforce: float | None,
    cnt_beadtype: str | None,
    cnt_functype: str | None,
    cnt_func_begin: int | None,
    cnt_func_end: int | None,
    cnt_base36: bool,
) -> tuple[Path, Path, list[Path]]:
    cnt_root = repo_cnt_root()
    script = cnt_root / "martini-cnt-generator.py"
    if not script.exists():
        raise FileNotFoundError(f"CNT generator not found: {script}")

    presets = {
        "cnt-m2": {
            "numrings": 12,
            "ringsize": 8,
            "bondlength": 0.47,
            "bondforce": 5000.0,
            "angleforce": 350.0,
            "beadtype": "CNP",
            "functype": "SNda",
            "func_begin": 1,
            "func_end": 1,
        },
        "cnt-m3": {
            "numrings": 12,
            "ringsize": 9,
            "bondlength": 0.41,
            "bondforce": 5000.0,
            "angleforce": 350.0,
            "beadtype": "SC5",
            "functype": "SNda",
            "func_begin": 1,
            "func_end": 1,
        },
    }
    preset = presets[mode]
    numrings = cnt_numrings if cnt_numrings is not None else preset["numrings"]
    ringsize = cnt_ringsize if cnt_ringsize is not None else preset["ringsize"]
    bondlength = cnt_bondlength if cnt_bondlength is not None else preset["bondlength"]
    bondforce = cnt_bondforce if cnt_bondforce is not None else preset["bondforce"]
    angleforce = cnt_angleforce if cnt_angleforce is not None else preset["angleforce"]
    beadtype = cnt_beadtype if cnt_beadtype is not None else preset["beadtype"]
    functype = cnt_functype if cnt_functype is not None else preset["functype"]
    func_begin = cnt_func_begin if cnt_func_begin is not None else preset["func_begin"]
    func_end = cnt_func_end if cnt_func_end is not None else preset["func_end"]

    gro_path = outdir / f"{basename}.gro"
    itp_path = outdir / f"{basename}.itp"
    posres_path = outdir / f"{basename}-posres.itp"
    cmd = [
        sys.executable,
        str(script),
        "-fn", str(outdir / basename),
        "-nr", str(numrings),
        "-rs", str(ringsize),
        "-bl", str(bondlength),
        "-bf", str(bondforce),
        "-af", str(angleforce),
        "-bt", beadtype,
        "-ft", functype,
        "-fb", str(func_begin),
        "-fe", str(func_end),
    ]
    if cnt_base36:
        cmd.append("--base36")

    result = subprocess.run(
        cmd,
        cwd=str(cnt_root),
        text=True,
        capture_output=True,
        check=False,
    )
    if result.returncode != 0:
        detail = (result.stderr or result.stdout or "unknown error").strip()
        raise RuntimeError(f"{mode} generator failed: {detail}")

    if not gro_path.exists():
        raise FileNotFoundError(f"Expected CNT structure not found: {gro_path}")
    if not itp_path.exists():
        raise FileNotFoundError(f"Expected CNT topology not found: {itp_path}")
    if not posres_path.exists():
        raise FileNotFoundError(f"Expected CNT position restraints not found: {posres_path}")

    rewrite_cnt_posres_include(itp_path, posres_path.name)
    return gro_path, itp_path, [posres_path]


def run_graphene_generator(
    mode: str,
    outdir: Path,
    basename: str,
    lx: float,
    ly: float,
    lz: float,
    graphite_layers: int,
    graphite_spacing: float,
) -> tuple[Path, Path, list[Path]]:
    graphene_root = repo_graphene_root()
    outdir.mkdir(parents=True, exist_ok=True)

    if mode == "graphene":
        script = graphene_root / "Periodic" / "martini3-graphene-periodic.py"
        cmd = [
            sys.executable,
            str(script),
            "-x", str(lx),
            "-y", str(ly),
            "-z", str(lz),
            "-o", basename,
            "--output-dir", str(outdir),
            "--force",
            "--quiet",
        ]
        gro_path = outdir / f"{basename}.gro"
        itp_path = outdir / f"{basename}.itp"
        support_files: list[Path] = []
    elif mode == "graphene-finite":
        script = graphene_root / "Non-Periodic" / "martini3-graphene-topology.py"
        cmd = [
            sys.executable,
            str(script),
            "-x", str(lx),
            "-y", str(ly),
            "-z", str(lz),
            "-o", basename,
            "--output-dir", str(outdir),
            "--force",
            "--quiet",
        ]
        gro_path = outdir / f"{basename}.gro"
        itp_path = outdir / f"{basename}.itp"
        support_files = []
    else:
        script = graphene_root / "graphite" / "build_graphite_slab.py"
        graphene_prefix = f"{basename}_layer"
        cmd = [
            sys.executable,
            str(script),
            "-x", str(lx),
            "-y", str(ly),
            "--layers", str(graphite_layers),
            "--spacing", str(graphite_spacing),
            "-o", basename,
            "--graphene-prefix", graphene_prefix,
            "--output-dir", str(outdir),
            "--force",
            "--quiet",
        ]
        gro_path = outdir / f"{basename}.gro"
        layer_itp = outdir / f"{graphene_prefix}.itp"
        itp_path = outdir / f"{basename}.itp"
        support_files = [outdir / "posre_GRA.itp"]

    result = subprocess.run(
        cmd,
        cwd=str(graphene_root),
        text=True,
        capture_output=True,
        check=False,
    )
    if result.returncode != 0:
        detail = (result.stderr or result.stdout or "unknown error").strip()
        raise RuntimeError(f"{mode} generator failed: {detail}")

    if mode == "graphite":
        if not layer_itp.exists():
            raise FileNotFoundError(f"Expected graphite topology not found: {layer_itp}")
        shutil.copy(layer_itp, itp_path)

    if not gro_path.exists():
        raise FileNotFoundError(f"Expected surface structure not found: {gro_path}")
    if not itp_path.exists():
        raise FileNotFoundError(f"Expected surface topology not found: {itp_path}")

    return gro_path, itp_path, support_files


def main(argv: Iterable[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Martini surface generator for generic lattices, graphene, and graphite."
    )

    # Base parameters
    parser.add_argument(
        "--mode",
        type=str,
        default="2-1",
        help="Surface mode: 2-1, 4-1, graphene, graphene-periodic, graphene-finite, graphite, cnt-m2, or cnt-m3.",
    )
    parser.add_argument(
        "--bead",
        nargs="+",
        default=["P4"],
        help="Martini bead type(s). Provide multiple values to cycle bead types by layer, e.g. --bead P4 C1.",
    )
    parser.add_argument("--dx", type=float, default=0.47, help="Bead spacing (2-1) or C-C distance (4-1)")
    parser.add_argument("--lx", type=float, required=True, help="Surface length in X (nm)")
    parser.add_argument("--ly", type=float, required=True, help="Surface length in Y (nm)")
    parser.add_argument("--lz", type=float, default=10.0, help="Box height in Z (nm)")
    parser.add_argument("--resname", type=str, default="SRF", help="Residue name")
    parser.add_argument("--output", type=str, default="surface", help="Output basename")
    parser.add_argument("--charge", type=float, default=0.0, help="Charge per bead")

    # Local hexagonal surface parameters
    parser.add_argument("--layers", type=int, default=1, help="Number of layers for local 2-1 / 4-1 surfaces.")
    parser.add_argument("--dist-z", type=float, default=0.382, help="Interlayer spacing in nm for local multilayer surfaces.")
    parser.add_argument("--graphite-layers", type=int, default=5, help="Number of layers in graphite mode")
    parser.add_argument("--graphite-spacing", type=float, default=0.382, help="Interlayer spacing in graphite mode (nm)")
    parser.add_argument("--cnt-numrings", type=int, help="Number of CNT rings.")
    parser.add_argument("--cnt-ringsize", type=int, help="Number of beads per CNT ring.")
    parser.add_argument("--cnt-bondlength", type=float, help="CNT bond length in nm.")
    parser.add_argument("--cnt-bondforce", type=float, help="CNT bond force constant.")
    parser.add_argument("--cnt-angleforce", type=float, help="CNT angle force constant.")
    parser.add_argument("--cnt-beadtype", type=str, help="CNT regular bead type.")
    parser.add_argument("--cnt-functype", type=str, help="CNT functionalized-end bead type.")
    parser.add_argument("--cnt-func-begin", type=int, help="Number of functionalized CNT rings at the beginning.")
    parser.add_argument("--cnt-func-end", type=int, help="Number of functionalized CNT rings at the end.")
    parser.add_argument("--cnt-base36", action="store_true", help="Use base36 atom naming for CNT mode.")

    args = parser.parse_args(list(argv) if argv is not None else None)
    mode = normalize_surface_mode(args.mode)
    if not mode:
        parser.error(f"Unsupported --mode '{args.mode}'. Valid modes: {', '.join(VALID_SURFACE_MODES)}")
    if args.lx <= 0 or args.ly <= 0 or args.lz <= 0:
        parser.error("--lx, --ly, and --lz must be positive.")
    if args.dx <= 0:
        parser.error("--dx must be positive.")
    if args.layers <= 0:
        parser.error("--layers must be positive.")
    if args.dist_z <= 0:
        parser.error("--dist-z must be positive.")
    if args.graphite_layers <= 0:
        parser.error("--graphite-layers must be positive.")
    if args.graphite_spacing <= 0:
        parser.error("--graphite-spacing must be positive.")
    if args.cnt_numrings is not None and args.cnt_numrings <= 0:
        parser.error("--cnt-numrings must be positive.")
    if args.cnt_ringsize is not None and args.cnt_ringsize <= 0:
        parser.error("--cnt-ringsize must be positive.")
    if args.cnt_bondlength is not None and args.cnt_bondlength <= 0:
        parser.error("--cnt-bondlength must be positive.")
    if args.cnt_bondforce is not None and args.cnt_bondforce <= 0:
        parser.error("--cnt-bondforce must be positive.")
    if args.cnt_angleforce is not None and args.cnt_angleforce <= 0:
        parser.error("--cnt-angleforce must be positive.")
    if args.cnt_func_begin is not None and args.cnt_func_begin < 0:
        parser.error("--cnt-func-begin must be >= 0.")
    if args.cnt_func_end is not None and args.cnt_func_end < 0:
        parser.error("--cnt-func-end must be >= 0.")

    outdir = os.path.dirname(args.output) or "."
    basename = os.path.basename(args.output)
    outdir_path = Path(outdir).expanduser().resolve()
    os.makedirs(outdir_path, exist_ok=True)

    layer_beads = normalize_layer_beads(args.bead)
    atoms: List[Tuple[float, float, float, str]] = []
    final_lx, final_ly = args.lx, args.ly

    # =========================================================
    # MODE 2-1: STANDARD HEXAGONAL MAPPING
    # =========================================================
    if mode == "2-1":
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

        for layer in range(args.layers):
            # Graphite-like ABAB stacking in the coarse-grained 2-1 lattice:
            # place odd layers over triangular hollow sites instead of directly
            # on top of occupied sites.
            shift_x, shift_y = graphite_ab_shift(layer, 0.5 * a, (math.sqrt(3) / 6) * a)
            z_pos = 3.0 + (layer * args.dist_z)
            bead = bead_for_layer(layer_beads, layer)

            for i in range(nx):
                for j in range(ny):
                    for x, y, _ in atoms_unit:
                        atoms.append((x + i * lx_cell + shift_x, y + j * ly_cell + shift_y, z_pos, bead))

    # =========================================================
    # MODE 4-1: HONEYCOMB CARBON MAPPING
    # =========================================================
    elif mode == "4-1":
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
            # Graphite-like ABAB (Bernal) stacking.
            shift_x, shift_y = graphite_ab_shift(layer, 0.0, d_cc)
            z_pos = 3.0 + (layer * args.dist_z)
            bead = bead_for_layer(layer_beads, layer)
            
            for i in range(nx):
                for j in range(ny):
                    for ax, ay, az in base_unit:
                        atoms.append((ax + i * lx_cell + shift_x, ay + j * ly_cell + shift_y, az + z_pos, bead))
    elif mode in {"cnt-m2", "cnt-m3"}:
        gro_path, itp_path, support_files = run_cnt_generator(
            mode=mode,
            outdir=outdir_path,
            basename=basename,
            cnt_numrings=args.cnt_numrings,
            cnt_ringsize=args.cnt_ringsize,
            cnt_bondlength=args.cnt_bondlength,
            cnt_bondforce=args.cnt_bondforce,
            cnt_angleforce=args.cnt_angleforce,
            cnt_beadtype=args.cnt_beadtype,
            cnt_functype=args.cnt_functype,
            cnt_func_begin=args.cnt_func_begin,
            cnt_func_end=args.cnt_func_end,
            cnt_base36=args.cnt_base36,
        )

        active_dst = resolve_activeitp_destination(outdir)
        active_dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(itp_path, active_dst)
        copy_support_files(support_files, active_dst.parent)
        print(f"✔ {mode} CNT surface generated from cnt-martini.")
        return
    else:
        gro_path, itp_path, support_files = run_graphene_generator(
            mode=mode,
            outdir=outdir_path,
            basename=basename,
            lx=args.lx,
            ly=args.ly,
            lz=args.lz,
            graphite_layers=args.graphite_layers,
            graphite_spacing=args.graphite_spacing,
        )

        active_dst = resolve_activeitp_destination(outdir)
        active_dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(itp_path, active_dst)
        copy_support_files(support_files, active_dst.parent)
        print(f"✔ {mode} surface generated from Martini3-Graphene.")
        return

    # ---------------------------------------------------------
    # File Writing
    # ---------------------------------------------------------
    gro_path = Path(outdir_path, f"{basename}.gro")
    with gro_path.open("w") as fgro:
        fgro.write(f"Surface {mode} | {args.lx}x{args.ly} nm | Layers: {args.layers}\n")
        fgro.write(f"{len(atoms):5d}\n")
        for i, (x, y, z, bead) in enumerate(atoms, start=1):
            fgro.write(f"{1:5d}{args.resname:<4}{bead:>4}{i:6d}{x:8.3f}{y:8.3f}{z:8.3f}\n")
        fgro.write(f"{final_lx:10.5f}{final_ly:10.5f}{args.lz:10.5f}\n")

    itp_path = Path(outdir_path, f"{basename}.itp")
    write_local_surface_itp(
        itp_path=itp_path,
        resname=args.resname,
        atoms=atoms,
        charge=args.charge,
    )

    # Backward-compatible copy for standalone tooling/tests.
    active_dst = resolve_activeitp_destination(outdir)
    active_dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(itp_path, active_dst)

    print(f"✔ {mode} Surface ({'Honeycomb Carbon' if mode == '4-1' else 'Hex Mapping'}) generated with {len(atoms)} beads.")

def resolve_activeitp_destination(outdir: str | Path) -> Path:
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
