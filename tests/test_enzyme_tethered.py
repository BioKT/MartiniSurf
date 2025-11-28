import numpy as np
from pathlib import Path
import surfmartini.enzyme_tethered as enz


# --------------------------------------------------------------
# Helper to create a minimal PDB
# --------------------------------------------------------------
def write_minimal_pdb(path: Path):
    pdb = """\
ATOM      1  CA  ALA A   1      10.000  10.000  10.000  1.00 20.00           C
ATOM      2  CA  ALA A   2      20.000  20.000  20.000  1.00 20.00           C
END
"""
    path.write_text(pdb)


# --------------------------------------------------------------
# Helper to create a minimal GRO for surface
# --------------------------------------------------------------
def write_minimal_gro(path: Path):
    gro = """\
Test GRO
    2
    1SURF   C1    1   0.000   0.000   0.000
    1SURF   C1    2   1.000   1.000   0.000
   2.00000   2.00000   2.00000
"""
    path.write_text(gro)


# --------------------------------------------------------------
# Test convert_pdb_to_gro
# --------------------------------------------------------------
def test_convert_pdb_to_gro(tmp_path):

    pdb = tmp_path / "test.pdb"
    gro = tmp_path / "test.gro"
    write_minimal_pdb(pdb)

    enz.convert_pdb_to_gro(pdb, gro)
    assert gro.exists()

    text = gro.read_text()
    assert "Converted from PDB" in text
    assert "2" in text.splitlines()[1]  # number of atoms


# --------------------------------------------------------------
# Test load_gro_coords
# --------------------------------------------------------------
def test_load_gro_coords(tmp_path):

    gro = tmp_path / "surf.gro"
    write_minimal_gro(gro)

    coords, atoms = enz.load_gro_coords(gro)

    assert coords.shape == (2, 3)
    assert len(atoms) == 2
    assert atoms[0][0] == 1   # resid
    assert atoms[0][1] == "SURF"


# --------------------------------------------------------------
# Test summarize_selected_residues
# --------------------------------------------------------------
def test_summarize_selected_residues(tmp_path):

    gro = tmp_path / "surf.gro"
    write_minimal_gro(gro)

    coords, atoms = enz.load_gro_coords(gro)

    # Residue 1 exists, residue 99 does not
    centroids = enz.summarize_selected_residues([1, 99], atoms, coords)

    assert centroids.shape == (1, 3)   # only 1 found
    assert np.allclose(centroids[0], coords.mean(axis=0))


# --------------------------------------------------------------
# Test auto_orient_from_anchor_residues
# --------------------------------------------------------------
def test_auto_orient_from_anchor_residues():

    enzyme = np.array([
        [0.0, 0.0, 5.0],
        [1.0, 1.0, 5.0],
    ])

    anchors = np.array([[0.5, 0.5, 5.0]])

    surface = np.array([
        [0.0, 0.0, 0.0],
        [2.0, 2.0, 0.0],
    ])

    result = enz.auto_orient_from_anchor_residues(
        enzyme, anchors, surface, target_z=5.0
    )

    assert result.shape == (2, 3)
    assert result[:, 2].min() >= 5.0   # enzyme must be above surface


# --------------------------------------------------------------
# Test save_full_system
# --------------------------------------------------------------
def test_save_full_system(tmp_path):

    out = tmp_path / "system.gro"

    # Minimal fake data
    surf_atoms = [(1, "SURF", "C1", 1), (1, "SURF", "C1", 2)]
    surf_coords = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 0.0]])

    enz_atoms = [(1, "ALA", "CA", 1), (2, "ALA", "CA", 2)]
    enz_coords = np.array([[10.0, 10.0, 10.0], [20.0, 20.0, 20.0]])

    enz.save_full_system(out, surf_atoms, surf_coords, enz_atoms, enz_coords)

    assert out.exists()
    text = out.read_text()
    assert "MartiniSurf oriented system" in text
    assert "4" in text.splitlines()[1]  # total atoms


