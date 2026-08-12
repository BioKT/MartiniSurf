from __future__ import annotations

import shlex
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class BuildConfig:
    input_path: Path | str = ""
    workdir: Path = Path(".")
    outdir: str = "Simulation_Files"
    moltype: str = "Protein"
    ff: str = "martini3001"
    dssp: bool = True
    go: bool = True
    go_eps: float = 9.415
    elastic: bool = False
    position_restraints: str = "backbone"
    pf: float = 1000.0
    maxwarn: int = 1
    merge_groups: list[str] = field(default_factory=list)
    surface_workflow: str = "Hexagonal Lattice"
    surface_mode: str = "4-1"
    surface_geometry: str = "planar"
    surface_path: Path | None = None
    lx: float = 6.0
    ly: float = 6.0
    dx: float = 0.47
    surface_beads: list[str] = field(default_factory=lambda: ["P4", "P4"])
    charge: float = 0.0
    surface_layers: int | None = 2
    surface_stacking: str = "hcp"
    surface_dist_z: float | None = None
    graphite_layers: int | None = None
    cnt_numrings: int | None = None
    cnt_ringsize: int | None = None
    orientation_mode: str = "Anchor"
    anchors: list[str] = field(default_factory=list)
    dist: float = 1.0
    ads_mode: bool = False
    balance_low_z: bool = False
    balance_low_z_fraction: float | None = None
    histag: bool = False
    linker_path: Path | None = None
    linker_groups: list[str] = field(default_factory=list)
    linker_prot_dist: float | None = None
    linker_surf_dist: float | None = None
    linker_protein_bead: str | None = None
    linker_surface_bead: str | None = None
    invert_linker: bool = False
    surface_linkers: int = 0
    substrate_path: Path | None = None
    substrate_itp_path: Path | None = None
    substrate_count: int = 0
    solvate: bool = True
    ionize: bool = True
    salt_conc: float = 0.15
    water_mix: str = ""
    water_mix_seed: int = 42
    martinize_extra_args: list[str] = field(default_factory=list)


def _append_value(args: list[str], flag: str, value: object | None) -> None:
    if value is None:
        return
    if isinstance(value, str) and not value.strip():
        return
    args.extend([flag, str(value)])


def _append_group(args: list[str], flag: str, value: str) -> None:
    tokens = shlex.split(value)
    if tokens:
        args.append(flag)
        args.extend(tokens)


def build_args(config: BuildConfig) -> list[str]:
    args = [
        "--pdb",
        str(config.input_path),
        "--moltype",
        config.moltype,
        "--ff",
        config.ff,
        "--p",
        config.position_restraints,
        "--pf",
        f"{config.pf:g}",
        "--maxwarn",
        str(config.maxwarn),
    ]

    args.extend([
        "--surface-mode",
        config.surface_mode,
        "--surface-geometry",
        config.surface_geometry,
        "--dx",
        f"{config.dx:g}",
        "--charge",
        f"{config.charge:g}",
        "--dist",
        f"{config.dist:g}",
        "--outdir",
        config.outdir,
    ])

    args.append("--dssp" if config.dssp else "--no-dssp")

    if config.go:
        args.append("--go")
        _append_value(args, "--go-eps", config.go_eps)
    if config.elastic:
        args.append("--elastic")

    if config.surface_path:
        _append_value(args, "--surface", config.surface_path)
    else:
        _append_value(args, "--lx", config.lx)
        _append_value(args, "--ly", config.ly)
        if config.surface_beads:
            args.append("--surface-bead")
            args.extend([bead for bead in config.surface_beads if bead])

    _append_value(args, "--surface-layers", config.surface_layers)
    _append_value(args, "--surface-stacking", config.surface_stacking)
    _append_value(args, "--surface-dist-z", config.surface_dist_z)
    _append_value(args, "--graphite-layers", config.graphite_layers)
    _append_value(args, "--cnt-numrings", config.cnt_numrings)
    _append_value(args, "--cnt-ringsize", config.cnt_ringsize)

    for merge_group in config.merge_groups:
        _append_value(args, "--merge", merge_group)

    needs_linker_file = config.orientation_mode == "Linker" or bool(config.surface_linkers)
    if needs_linker_file:
        _append_value(args, "--linker", config.linker_path)
        _append_value(args, "--linker-surface-bead", config.linker_surface_bead)
    if config.orientation_mode == "Linker":
        for group in config.linker_groups:
            _append_group(args, "--linker-group", group)
        _append_value(args, "--linker-prot-dist", config.linker_prot_dist)
        _append_value(args, "--linker-surf-dist", config.linker_surf_dist)
        _append_value(args, "--linker-protein-bead", config.linker_protein_bead)
        if config.invert_linker:
            args.append("--invert-linker")
    else:
        for anchor in config.anchors:
            _append_group(args, "--anchor", anchor)
        if config.ads_mode:
            args.append("--ads-mode")
        if config.balance_low_z:
            args.append("--balance-low-z")
            _append_value(args, "--balance-low-z-fraction", config.balance_low_z_fraction)
        if config.histag:
            args.append("--histag")
    if config.surface_linkers:
        _append_value(args, "--surface-linkers", config.surface_linkers)

    if config.substrate_path and config.substrate_count > 0:
        _append_value(args, "--substrate", config.substrate_path)
        _append_value(args, "--substrate-itp", config.substrate_itp_path)
        _append_value(args, "--substrate-count", config.substrate_count)

    if config.solvate:
        args.append("--solvate")
    if config.ionize:
        args.append("--ionize")
        _append_value(args, "--salt-conc", config.salt_conc)
    if config.water_mix.strip():
        _append_value(args, "--water-mix", config.water_mix.strip())
        _append_value(args, "--water-mix-seed", config.water_mix_seed)

    for extra in config.martinize_extra_args:
        if extra.strip():
            args.extend(["--martinize-extra-args", extra.strip()])

    return args


def shell_command(args: list[str]) -> str:
    return "python -m martinisurf " + " ".join(shlex.quote(arg) for arg in args)
