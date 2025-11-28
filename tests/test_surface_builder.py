import os
import math
from surfmartini.surface_builder import main


# ============================================================
# BASIC BUILD TEST
# ============================================================
def test_surface_builder_basic(tmp_path):
    out = tmp_path / "mysurface"

    main([
        "--lx", "5",
        "--ly", "5",
        "--dx", "0.5",
        "--output", str(out)
    ])

    gro = tmp_path / "mysurface.gro"
    itp = tmp_path / "mysurface.itp"
    dst = tmp_path / "ActiveITP" / "surface.itp"

    assert gro.exists()
    assert itp.exists()
    assert dst.exists()

    text = gro.read_text()
    lines = text.strip().splitlines()

    # Check title
    assert "surface" in lines[0].lower()

    # Atom count
    n_atoms = int(lines[1])
    assert n_atoms > 0

    # Box vector line
    box_line = lines[-1].split()
    assert len(box_line) == 3


# ============================================================
# PARAMETERS TEST
# ============================================================
def test_surface_parameters(tmp_path):
    out = tmp_path / "surfX"

    main([
        "--lx", "10",
        "--ly", "7",
        "--dx", "0.4",
        "--bead", "Qd",
        "--resname", "ABC",
        "--charge", "0.55",
        "--output", str(out)
    ])

    itp_file = tmp_path / "surfX.itp"
    gro_file = tmp_path / "surfX.gro"

    itp = itp_file.read_text()
    gro = gro_file.read_text()

    assert "ABC" in gro
    assert "Qd" in gro

    assert "0.550" in itp
    assert "Qd" in itp


# ============================================================
# LATTICE ATOM COUNT TEST
# ============================================================
def test_lattice_atom_count(tmp_path):
    out = tmp_path / "surf"

    dx = 0.47
    lx = 6.0
    ly = 6.0

    main([
        "--lx", str(lx),
        "--ly", str(ly),
        "--dx", str(dx),
        "--output", str(out)
    ])

    gro_file = tmp_path / "surf.gro"
    text = gro_file.read_text()
    n_atoms = int(text.splitlines()[1])

    # Recompute lattice parameters exactly like surface_builder
    scale = dx / 0.142
    a = 0.246 * scale
    Lx_cell = 3.0 * a
    Ly_cell = math.sqrt(3) * a

    nx = max(1, round(lx / Lx_cell))
    ny = max(1, round(ly / Ly_cell))

    expected = nx * ny * 6  # 6 beads per hex cell

    assert n_atoms == expected


# ============================================================
# STANDALONE MODE
# ============================================================
def test_standalone_mode_copy(tmp_path):
    out = tmp_path / "mysurf"

    main([
        "--lx", "4",
        "--ly", "4",
        "--output", str(out)
    ])

    dst = tmp_path / "ActiveITP" / "surface.itp"
    assert dst.exists()

