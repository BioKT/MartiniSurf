import numpy as np
from pathlib import Path
import martinisurf.system_tethered as enz


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


def write_system_gro(path: Path):
    gro = """\
System GRO
    2
    1MOL    B1    1   0.500   0.500   1.000
    1MOL    B2    2   0.700   0.700   1.100
   2.00000   2.00000   2.00000
"""
    path.write_text(gro)


def write_linker_gro(path: Path):
    gro = """\
Linker GRO
    3
    1LNK    L1    1   0.000   0.000   0.000
    1LNK    L2    2   0.000   0.000   0.100
    1LNK    L3    3   0.000   0.000   0.200
   2.00000   2.00000   2.00000
"""
    path.write_text(gro)


def atom_z_by_name(gro_path: Path, atom_name: str) -> float:
    for line in gro_path.read_text().splitlines()[2:-1]:
        if line[10:15].strip() == atom_name:
            return float(line[36:44])
    raise ValueError(f"Atom {atom_name} not found")


def linker_atom_sequence(gro_path: Path) -> list[str]:
    names = []
    for line in gro_path.read_text().splitlines()[2:-1]:
        if line[5:10].strip() == "LNK":
            names.append(line[10:15].strip())
    return names


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


def test_auto_orient_prevents_surface_penetration_when_anchors_are_high():
    enzyme = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 1.0, 100.0],
    ])
    anchors = np.array([
        [1.0, 1.0, 100.0],
    ])
    surface = np.array([
        [0.0, 0.0, 0.0],
        [2.0, 2.0, 0.0],
    ])

    result = enz.auto_orient_from_anchor_residues(
        enzyme, anchors, surface, target_z=10.0
    )

    # Even with very high anchors, the full system must stay above the surface.
    assert result[:, 2].min() >= 1.0


def test_auto_orient_respects_anchor_distance_after_flip():
    # Geometry chosen so anchor normal can trigger the flip branch.
    enzyme = np.array([
        [0.0, 0.0, 0.0],   # anchor atom
        [1.0, 0.0, 0.0],   # anchor atom
        [0.0, 1.0, 5.0],   # high point
        [0.0, 0.0, -5.0],  # low point
    ])
    anchors = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
    ])
    surface = np.array([
        [0.0, 0.0, 0.0],
        [2.0, 2.0, 0.0],
    ])

    target = 10.0
    result = enz.auto_orient_from_anchor_residues(
        enzyme, anchors, surface, target_z=target
    )

    # Anchor atoms are indices 0 and 1 in this synthetic setup.
    transformed_anchor_min = result[[0, 1], 2].min()
    # Allow a small tolerance due to rotation/format precision.
    assert transformed_anchor_min >= target - 1e-6


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

    # NEW: explicit box line (same format as GRO last line)
    box_line = "2.00000 2.00000 2.00000"

    enz.save_full_system(
        out,
        surf_atoms,
        surf_coords,
        enz_atoms,
        enz_coords,
        box_line,
    )

    assert out.exists()
    text = out.read_text()
    lines = text.splitlines()

    assert "MartiniSurf oriented system" in lines[0]
    assert int(lines[1]) == 4  # total atoms
    box_vals = lines[-1].split()
    assert len(box_vals) >= 3
    # XY preserved from input box line.
    assert abs(float(box_vals[0]) - 2.0) < 1e-6
    assert abs(float(box_vals[1]) - 2.0) < 1e-6
    # Z must be at least max(enz_z) + 3 nm.
    assert float(box_vals[2]) >= 5.0


def test_save_full_system_expands_xy_box_when_needed(tmp_path):
    out = tmp_path / "expanded_xy.gro"

    surf_atoms = [(1, "SRF", "C1", 1)]
    surf_coords = np.array([[0.0, 0.0, 0.0]])
    enz_atoms = [(1, "ALA", "CA", 1)]
    enz_coords = np.array([[200.0, 210.0, 20.0]])  # 20.0 nm, 21.0 nm
    box_line = "2.00000 2.00000 2.00000"

    enz.save_full_system(
        out,
        surf_atoms,
        surf_coords,
        enz_atoms,
        enz_coords,
        box_line,
    )

    lines = out.read_text().splitlines()
    box_vals = lines[-1].split()
    assert float(box_vals[0]) >= 20.0
    assert float(box_vals[1]) >= 21.0


def test_invert_linker_switches_attachment_side(tmp_path):
    surface = tmp_path / "surface.gro"
    system = tmp_path / "system.gro"
    linker = tmp_path / "linker.gro"
    out_default = tmp_path / "out_default.gro"
    out_invert = tmp_path / "out_invert.gro"

    write_minimal_gro(surface)
    write_system_gro(system)
    write_linker_gro(linker)

    enz.main([
        "--surface", str(surface),
        "--system", str(system),
        "--out", str(out_default),
        "--linker-gro", str(linker),
        "--linker-group", "1", "1",
        "--linker-prot-dist", "0.0",
        "--linker-surf-dist", "0.0",
    ])

    enz.main([
        "--surface", str(surface),
        "--system", str(system),
        "--out", str(out_invert),
        "--linker-gro", str(linker),
        "--linker-group", "1", "1",
        "--linker-prot-dist", "0.0",
        "--linker-surf-dist", "0.0",
        "--invert-linker",
    ])

    assert linker_atom_sequence(out_default) == ["L1", "L2", "L3"]
    assert linker_atom_sequence(out_invert) == ["L3", "L2", "L1"]


def test_linker_mode_raises_clear_error_for_missing_group_residues(tmp_path):
    surface = tmp_path / "surface.gro"
    system = tmp_path / "system.gro"
    linker = tmp_path / "linker.gro"
    out = tmp_path / "out_bad_group.gro"

    write_minimal_gro(surface)
    write_system_gro(system)
    write_linker_gro(linker)

    try:
        enz.main([
            "--surface", str(surface),
            "--system", str(system),
            "--out", str(out),
            "--linker-gro", str(linker),
            "--linker-group", "1", "99",
            "--linker-prot-dist", "0.0",
            "--linker-surf-dist", "0.0",
        ])
    except ValueError as exc:
        msg = str(exc)
        assert "No atoms found for linker groups" in msg
        assert "Requested residues: [99]" in msg
    else:
        raise AssertionError("Expected ValueError for missing linker-group residues.")

