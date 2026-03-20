from pathlib import Path

from martinisurf.pipeline import (
    _build_generated_surface_args,
    _effective_surface_geometry,
    _resolve_generated_surface_mode,
    build_parser,
)


def test_surface_mode_defaults_to_2_1():
    parser = build_parser()
    args = parser.parse_args([
        "--pdb", "1RJW",
        "--anchor", "1", "1",
        "--lx", "10",
        "--ly", "10",
    ])
    assert args.surface_mode == "2-1"


def test_surface_mode_accepts_4_1():
    parser = build_parser()
    args = parser.parse_args([
        "--pdb", "1RJW",
        "--anchor", "1", "1",
        "--lx", "10",
        "--ly", "10",
        "--surface-mode", "4-1",
    ])
    assert args.surface_mode == "4-1"


def test_surface_mode_accepts_graphene():
    parser = build_parser()
    args = parser.parse_args([
        "--pdb", "1RJW",
        "--anchor", "1", "1",
        "--lx", "10",
        "--ly", "10",
        "--surface-mode", "graphene",
    ])
    assert args.surface_mode == "graphene"


def test_surface_mode_accepts_cnt_alias():
    parser = build_parser()
    args = parser.parse_args([
        "--pdb", "1RJW",
        "--anchor", "1", "1",
        "--lx", "10",
        "--ly", "10",
        "--surface-mode", "cnt",
        "--cnt-numrings", "24",
        "--cnt-ringsize", "12",
    ])
    assert args.surface_mode == "cnt"
    assert args.cnt_numrings == 24
    assert args.cnt_ringsize == 12


def test_surface_geometry_defaults_to_planar():
    parser = build_parser()
    args = parser.parse_args([
        "--pdb", "1RJW",
        "--anchor", "1", "1",
        "--lx", "10",
        "--ly", "10",
    ])
    assert args.surface_geometry == "planar"


def test_surface_geometry_accepts_3d():
    parser = build_parser()
    args = parser.parse_args([
        "--pdb", "1RJW",
        "--anchor", "1", "1",
        "--lx", "10",
        "--ly", "10",
        "--surface-geometry", "3d",
    ])
    assert args.surface_geometry == "3d"


def test_surface_mode_accepts_graphite_parameters():
    parser = build_parser()
    args = parser.parse_args([
        "--pdb", "1RJW",
        "--anchor", "1", "1",
        "--lx", "10",
        "--ly", "10",
        "--surface-mode", "graphite",
        "--graphite-layers", "4",
        "--graphite-spacing", "0.382",
    ])
    assert args.surface_mode == "graphite"
    assert args.graphite_layers == 4
    assert args.graphite_spacing == 0.382


def test_cnt_alias_resolves_to_martini3_for_protein():
    parser = build_parser()
    args = parser.parse_args([
        "--pdb", "1RJW",
        "--anchor", "1", "1",
        "--lx", "10",
        "--ly", "10",
        "--surface-mode", "cnt",
    ])
    assert _resolve_generated_surface_mode(args) == "cnt-m3"


def test_cnt_alias_resolves_to_martini2_for_dna():
    parser = build_parser()
    args = parser.parse_args([
        "--dna",
        "--pdb", "4C64",
        "--anchor", "1", "1",
        "--lx", "10",
        "--ly", "10",
        "--surface-mode", "cnt",
    ])
    assert _resolve_generated_surface_mode(args) == "cnt-m2"


def test_cnt_surface_uses_3d_geometry_automatically():
    parser = build_parser()
    args = parser.parse_args([
        "--pdb", "1RJW",
        "--anchor", "1", "1",
        "--lx", "10",
        "--ly", "10",
        "--surface-mode", "cnt",
    ])
    assert _effective_surface_geometry(args) == "3d"


def test_external_surface_uses_3d_geometry_automatically():
    parser = build_parser()
    args = parser.parse_args([
        "--pdb", "1RJW",
        "--anchor", "1", "1",
        "--surface", "custom_surface.gro",
    ])
    assert _effective_surface_geometry(args) == "3d"


def test_external_surface_forces_3d_even_if_planar_requested():
    parser = build_parser()
    args = parser.parse_args([
        "--pdb", "1RJW",
        "--anchor", "1", "1",
        "--surface", "custom_surface.gro",
        "--surface-geometry", "planar",
    ])
    assert _effective_surface_geometry(args) == "3d"


def test_generated_surface_args_forward_cnt_flags():
    parser = build_parser()
    args = parser.parse_args([
        "--dna",
        "--pdb", "4C64",
        "--anchor", "1", "1",
        "--lx", "10",
        "--ly", "10",
        "--surface-mode", "cnt",
        "--cnt-numrings", "20",
        "--cnt-ringsize", "11",
        "--cnt-bondlength", "0.5",
        "--cnt-bondforce", "6000",
        "--cnt-angleforce", "400",
        "--cnt-beadtype", "CNP",
        "--cnt-functype", "SNda",
        "--cnt-func-begin", "2",
        "--cnt-func-end", "3",
        "--cnt-base36",
    ])
    builder_args = _build_generated_surface_args(args, Path("surface"))
    assert builder_args[:4] == ["--mode", "cnt-m2", "--lx", "10.0"]
    assert "--martini-version" in builder_args
    mv_index = builder_args.index("--martini-version")
    assert builder_args[mv_index + 1] == "2"
    assert "--cnt-numrings" in builder_args
    assert "--cnt-ringsize" in builder_args
    assert "--cnt-bondlength" in builder_args
    assert "--cnt-bondforce" in builder_args
    assert "--cnt-angleforce" in builder_args
    assert "--cnt-beadtype" in builder_args
    assert "--cnt-functype" in builder_args
    assert "--cnt-func-begin" in builder_args
    assert "--cnt-func-end" in builder_args
    assert "--cnt-base36" in builder_args


def test_surface_bead_accepts_multiple_values_and_forwards_them():
    parser = build_parser()
    args = parser.parse_args([
        "--pdb", "1RJW",
        "--anchor", "1", "1",
        "--lx", "10",
        "--ly", "10",
        "--surface-mode", "2-1",
        "--surface-layers", "2",
        "--surface-bead", "P4", "C1",
    ])
    assert args.surface_bead == ["P4", "C1"]

    builder_args = _build_generated_surface_args(args, Path("surface"))
    bead_index = builder_args.index("--bead")
    mv_index = builder_args.index("--martini-version")
    assert builder_args[mv_index + 1] == "3"
    assert builder_args[bead_index + 1:bead_index + 3] == ["P4", "C1"]
