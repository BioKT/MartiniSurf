from pathlib import Path

from streamlit_app.command_builder import BuildConfig, build_args
from streamlit_app.linker_generator import parse_linker_gro, parse_martini_mapper_report, safe_molecule_name


def test_parse_linker_gro_reads_bead_labels():
    gro_path = Path("martinisurf/examples/protein/03_linker_surface_decoration/inputs/EPOXY.gro")

    beads = parse_linker_gro(gro_path)

    assert [bead.label for bead in beads] == ["1: P01", "2: P02"]
    assert beads[0].resname == "EPOX"


def test_parse_linker_gro_accepts_compact_martini_mapper_lines(tmp_path):
    gro_path = tmp_path / "LINKER.gro"
    gro_path.write_text(
        "1LINKER\n"
        "1\n"
        "    1res  C1    1   -0.000  -0.000  -0.000\n"
        "  10.00000  10.00000  10.00000\n"
    )

    beads = parse_linker_gro(gro_path)

    assert [bead.label for bead in beads] == ["1: C1"]
    assert beads[0].resname == "res"


def test_safe_molecule_name_keeps_generator_output_predictable():
    assert safe_molecule_name(" linker 01!") == "linker_01"
    assert safe_molecule_name("   ") == "LINKER"


def test_parse_martini_mapper_report_explains_bead_chemistry(tmp_path):
    gro_path = tmp_path / "LINKER.gro"
    gro_path.write_text(
        "1LINKER\n"
        "3\n"
        "    1res  C1    1    0.048  -0.169  -0.011\n"
        "    1res  C2    2   -0.208   0.047   0.009\n"
        "    1res  C3    3    0.118   0.131   0.004\n"
        "  10.00000  10.00000  10.00000\n"
    )
    report_path = tmp_path / "LINKER.txt"
    report_path.write_text(
        "Compound: LINKER\n"
        "SMILES:   Oc1ccccc1O\n\n"
        "All beads\n"
        "--------------------------------------------------------------------------------\n"
        "index    bead_type  atoms_idx            atoms_elem           #conns       mass     x        y        z\n"
        "0        SN6        0,1,2                O,c,c                2            54.0     0.0477   -0.1685  -0.0107\n"
        "1        TC5        3,4                  c,c                  2            36.0     -0.2077  0.0469   0.0088\n"
        "2        SN6        5,6,7                c,c,O                2            54.0     0.1185   0.1310   0.0037\n\n"
        "Bead connections\n"
    )

    rows = parse_martini_mapper_report(report_path, parse_linker_gro(gro_path))

    assert rows[0].bead_label == "1: C1"
    assert rows[0].martini_type == "SN6"
    assert rows[0].atom_indices_1based == [1, 2, 3]
    assert rows[2].atom_indices_1based == [6, 7, 8]


def test_build_args_passes_selected_linker_beads():
    args = build_args(
        BuildConfig(
            input_path="1UBQ",
            orientation_mode="Linker",
            linker_path=Path("linker.gro"),
            linker_groups=["A 73"],
            linker_protein_bead="1: C1",
            linker_surface_bead="2: C2",
        )
    )

    assert "--linker-protein-bead" in args
    assert args[args.index("--linker-protein-bead") + 1] == "1: C1"
    assert "--linker-surface-bead" in args
    assert args[args.index("--linker-surface-bead") + 1] == "2: C2"
