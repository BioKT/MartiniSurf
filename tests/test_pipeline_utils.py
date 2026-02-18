from pathlib import Path

from martinisurf.pipeline import _backup_existing_output_dir, _write_minimal_surface_itp


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
