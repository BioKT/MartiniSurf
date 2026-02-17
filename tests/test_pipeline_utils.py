from pathlib import Path

from martinisurf.pipeline import _backup_existing_output_dir


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
