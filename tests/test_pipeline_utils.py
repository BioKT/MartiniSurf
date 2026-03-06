from pathlib import Path

import pytest

from martinisurf.pipeline import _backup_existing_output_dir, _write_minimal_surface_itp
from martinisurf import pipeline
import martinisurf


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
    assert cfg["anchor_landmark_mode"] == "residue"
    assert cfg["balance_low_z"] is False
    assert abs(cfg["balance_low_z_fraction"] - 0.2) < 1e-12
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
    assert cfg["anchor_landmark_mode"] == "group"


def test_load_pre_cg_complex_config_reads_low_z_balance_options(tmp_path):
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
        "  balance_low_z: true\n"
        "  balance_low_z_fraction: 0.35\n"
        "cofactor:\n"
        "  molname: COF\n"
        "  itp: cofactor.itp\n"
        "  count: 1\n"
        "topology:\n"
        "  protein_itp: protein.itp\n"
        "  include_go: false\n"
    )

    cfg = pipeline._load_pre_cg_complex_config(input_dir / "complex_config.yaml")
    assert cfg["balance_low_z"] is True
    assert abs(cfg["balance_low_z_fraction"] - 0.35) < 1e-12


def test_load_pre_cg_complex_config_rejects_invalid_low_z_fraction(tmp_path):
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
        "  balance_low_z: true\n"
        "  balance_low_z_fraction: 2.0\n"
        "cofactor:\n"
        "  molname: COF\n"
        "  itp: cofactor.itp\n"
        "  count: 1\n"
        "topology:\n"
        "  protein_itp: protein.itp\n"
        "  include_go: false\n"
    )

    with pytest.raises(ValueError, match="balance_low_z_fraction"):
        pipeline._load_pre_cg_complex_config(input_dir / "complex_config.yaml")


def test_load_pre_cg_complex_config_resolves_chain_based_anchor_groups_with_reference_pdb(tmp_path):
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    (input_dir / "complex.gro").write_text("X\n0\n1 1 1\n")
    (input_dir / "protein.itp").write_text("[ moleculetype ]\nPROT 1\n")
    (input_dir / "cofactor.itp").write_text("[ moleculetype ]\nCOF 1\n")
    (input_dir / "protein_ref.pdb").write_text(
        "ATOM      1  N   GLY A   8       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  N   SER A  10       0.100   0.000   0.000  1.00  0.00           N\n"
        "ATOM      3  N   THR A  11       0.200   0.000   0.000  1.00  0.00           N\n"
        "ATOM      4  N   GLY D   8       0.300   0.000   0.000  1.00  0.00           N\n"
        "ATOM      5  N   SER D  10       0.400   0.000   0.000  1.00  0.00           N\n"
        "ATOM      6  N   THR D  11       0.500   0.000   0.000  1.00  0.00           N\n"
    )
    (input_dir / "complex_config.yaml").write_text(
        "mode: pre_cg_complex\n"
        "complex_gro: complex.gro\n"
        "protein:\n"
        "  molname: PROT\n"
        "  reference_pdb: protein_ref.pdb\n"
        "  anchor_groups: [\"A 8 10 11\", \"D 8 10 11\"]\n"
        "cofactor:\n"
        "  molname: COF\n"
        "  itp: cofactor.itp\n"
        "  count: 1\n"
        "topology:\n"
        "  protein_itp: protein.itp\n"
        "  include_go: false\n"
    )

    cfg = pipeline._load_pre_cg_complex_config(input_dir / "complex_config.yaml")
    assert cfg["anchor_groups"] == [[1, 1, 2, 3], [2, 4, 5, 6]]
    assert cfg["anchor_landmark_mode"] == "group"
    assert cfg["reference_pdb"] == (input_dir / "protein_ref.pdb").resolve()


def test_load_pre_cg_complex_config_requires_reference_pdb_for_chain_based_anchor_groups(tmp_path):
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
        "  anchor_groups: [\"A 8 10 11\"]\n"
        "cofactor:\n"
        "  molname: COF\n"
        "  itp: cofactor.itp\n"
        "  count: 1\n"
        "topology:\n"
        "  protein_itp: protein.itp\n"
        "  include_go: false\n"
    )

    with pytest.raises(ValueError, match="protein.reference_pdb"):
        pipeline._load_pre_cg_complex_config(input_dir / "complex_config.yaml")


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
    assert cfg["anchor_landmark_mode"] == "residue"


def test_normalize_cli_residue_groups_resolves_chain_based_anchor_groups(tmp_path):
    pdb = tmp_path / "cleaned_input.pdb"
    pdb.write_text(
        "ATOM      1  N   GLY A   8       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  GLY A   8       0.100   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  N   SER A  10       0.200   0.000   0.000  1.00  0.00           N\n"
        "ATOM      4  N   THR D   8       0.300   0.000   0.000  1.00  0.00           N\n"
        "ATOM      5  N   TYR D  10       0.400   0.000   0.000  1.00  0.00           N\n"
    )

    groups = pipeline._normalize_cli_residue_groups(
        [["A", "8", "10"], ["D", "8", "10"]],
        pdb,
        "--anchor",
    )

    assert groups == [[1, 1, 2], [2, 3, 4]]


def test_normalize_cli_residue_groups_keeps_legacy_numeric_syntax(tmp_path):
    pdb = tmp_path / "cleaned_input.pdb"
    pdb.write_text(
        "ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00  0.00           N\n"
    )

    groups = pipeline._normalize_cli_residue_groups(
        [["7", "11", "13"]],
        pdb,
        "--anchor",
    )

    assert groups == [[7, 11, 13]]


def test_normalize_cli_residue_groups_raises_clear_error_for_missing_chain_residue(tmp_path):
    pdb = tmp_path / "cleaned_input.pdb"
    pdb.write_text(
        "ATOM      1  N   GLY A   8       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  N   THR D   8       0.300   0.000   0.000  1.00  0.00           N\n"
    )

    with pytest.raises(ValueError, match="could not resolve D 10"):
        pipeline._normalize_cli_residue_groups(
            [["D", "10"]],
            pdb,
            "--linker-group",
        )


def test_normalize_ion_atom_names_from_itp_converts_legacy_plus_minus_labels(tmp_path):
    top_dir = tmp_path / "0_topology"
    itp_dir = top_dir / "system_itp"
    itp_dir.mkdir(parents=True)

    src_ions = Path(martinisurf.__file__).resolve().parent / "system_itp" / "martini_v2.0_ions.itp"
    (itp_dir / "martini_v2.0_ions.itp").write_text(src_ions.read_text())

    gro = tmp_path / "final_system.gro"
    pipeline._write_gro_records(
        str(gro),
        "ions test",
        [
            {"resid": 1, "resname": "NA+", "atomname": "NA+", "atomid": 1, "x": 0.1, "y": 0.1, "z": 0.1},
            {"resid": 2, "resname": "CL-", "atomname": "CL-", "atomid": 2, "x": 0.2, "y": 0.2, "z": 0.2},
        ],
        [1.0, 1.0, 1.0],
    )

    changed = pipeline._normalize_ion_atom_names_from_itp(top_dir=top_dir, gro_path=gro)
    assert changed is True

    _, records, _ = pipeline._read_gro_records(str(gro))
    assert records[0]["resname"] == "NA"
    assert records[0]["atomname"] == "NA"
    assert records[1]["resname"] == "CL"
    assert records[1]["atomname"] == "CL"
