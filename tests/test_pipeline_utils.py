from pathlib import Path

import pytest

from martinisurf.pipeline import _backup_existing_output_dir, _write_minimal_surface_itp
from martinisurf import pipeline


def test_backup_existing_output_dir_moves_directory(tmp_path):
    outdir = tmp_path / "Simulation_Files"
    outdir.mkdir()
    marker = outdir / "keep.txt"
    marker.write_text("important")

    backup = _backup_existing_output_dir(outdir)

    assert backup is not None
    assert not outdir.exists()
    assert backup.exists()
    assert (backup / "keep.txt").read_text() == "important"


def test_backup_existing_output_dir_noop_when_missing(tmp_path):
    outdir = tmp_path / "Missing_Output"
    backup = _backup_existing_output_dir(outdir)
    assert backup is None


def test_write_minimal_surface_itp_writes_expected_blocks(tmp_path):
    itp = tmp_path / "surface.itp"
    _write_minimal_surface_itp(itp, resname="SRF", bead="C1", charge=0.0)
    text = itp.read_text()
    assert "[ moleculetype ]" in text
    assert "[ atoms ]" in text
    assert "SRF" in text
    assert "C1" in text


def test_load_pre_cg_complex_config_reads_required_fields(tmp_path):
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    (input_dir / "complex.gro").write_text("X\n0\n1 1 1\n")
    (input_dir / "protein.itp").write_text("[ moleculetype ]\nPROT 1\n")
    (input_dir / "cofactor.itp").write_text("[ moleculetype ]\nCOF 1\n")
    (input_dir / "go_atomtypes.itp").write_text("; go\n")
    (input_dir / "go_nbparams.itp").write_text("; go\n")
    (input_dir / "complex_config.yaml").write_text(
        "mode: pre_cg_complex\n"
        "complex_gro: complex.gro\n"
        "protein:\n"
        "  molname: PROT\n"
        "  orient_by_residues: [45, 67, 120]\n"
        "cofactor:\n"
        "  molname: COF\n"
        "  itp: cofactor.itp\n"
        "  count: 2\n"
        "topology:\n"
        "  protein_itp: protein.itp\n"
        "  include_go: true\n"
        "  go_files_glob: \"go_*\"\n"
    )

    cfg = pipeline._load_pre_cg_complex_config(input_dir / "complex_config.yaml")
    assert cfg["protein_molname"] == "PROT"
    assert cfg["cofactor_molname"] == "COF"
    assert cfg["cofactor_count"] == 2
    assert cfg["anchor_groups"] == [[1, 45, 67, 120]]
    assert len(cfg["go_files"]) == 2


def test_load_pre_cg_complex_config_requires_cofactor_count(tmp_path):
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    (input_dir / "complex.gro").write_text("X\n0\n1 1 1\n")
    (input_dir / "protein.itp").write_text("[ moleculetype ]\nPROT 1\n")
    (input_dir / "cofactor.itp").write_text("[ moleculetype ]\nCOF 1\n")
    (input_dir / "complex_config.yaml").write_text(
        "mode: pre_cg_complex\n"
        "complex_gro: complex.gro\n"
        "protein:\n"
        "  molname: PROT\n"
        "  orient_by_residues: [1]\n"
        "cofactor:\n"
        "  molname: COF\n"
        "  itp: cofactor.itp\n"
        "topology:\n"
        "  protein_itp: protein.itp\n"
        "  include_go: false\n"
    )

    with pytest.raises(ValueError, match="cofactor.count"):
        pipeline._load_pre_cg_complex_config(input_dir / "complex_config.yaml")


def test_load_pre_cg_complex_config_accepts_anchor_groups(tmp_path):
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    (input_dir / "complex.gro").write_text("X\n0\n1 1 1\n")
    (input_dir / "protein.itp").write_text("[ moleculetype ]\nPROT 1\n")
    (input_dir / "cofactor.itp").write_text("[ moleculetype ]\nCOF 1\n")
    (input_dir / "complex_config.yaml").write_text(
        "mode: pre_cg_complex\n"
        "complex_gro: complex.gro\n"
        "protein:\n"
        "  molname: PROT\n"
        "  anchor_groups: [\"1 8 10 11\", \"2 1025 1027 1028\"]\n"
        "cofactor:\n"
        "  molname: COF\n"
        "  itp: cofactor.itp\n"
        "  count: 1\n"
        "topology:\n"
        "  protein_itp: protein.itp\n"
        "  include_go: false\n"
    )

    cfg = pipeline._load_pre_cg_complex_config(input_dir / "complex_config.yaml")
    assert cfg["anchor_groups"] == [[1, 8, 10, 11], [2, 1025, 1027, 1028]]


def test_load_pre_cg_complex_config_allows_missing_anchor_groups(tmp_path):
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    (input_dir / "complex.gro").write_text("X\n0\n1 1 1\n")
    (input_dir / "protein.itp").write_text("[ moleculetype ]\nPROT 1\n")
    (input_dir / "cofactor.itp").write_text("[ moleculetype ]\nCOF 1\n")
    (input_dir / "complex_config.yaml").write_text(
        "mode: pre_cg_complex\n"
        "complex_gro: complex.gro\n"
        "protein:\n"
        "  molname: PROT\n"
        "cofactor:\n"
        "  molname: COF\n"
        "  itp: cofactor.itp\n"
        "  count: 1\n"
        "topology:\n"
        "  protein_itp: protein.itp\n"
        "  include_go: false\n"
    )

    cfg = pipeline._load_pre_cg_complex_config(input_dir / "complex_config.yaml")
    assert cfg["anchor_groups"] == []
