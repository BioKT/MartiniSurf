from pathlib import Path

from streamlit_app.linker_generator import parse_linker_gro, safe_molecule_name


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
