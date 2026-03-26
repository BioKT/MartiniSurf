import os
import math
from martinisurf.surface_builder import SIGMA_SPACING_FACTOR, main


def _itp_section_lines(text: str, section_name: str) -> list[str]:
    target = f"[ {section_name} ]"
    in_section = False
    lines: list[str] = []
    for raw in text.splitlines():
        line = raw.strip()
        if line.startswith("["):
            in_section = line.lower() == target
            continue
        if in_section and line and not line.startswith(";"):
            lines.append(line)
    return lines


def _bond_degree_counts(bond_lines: list[str]) -> dict[int, int]:
    counts: dict[int, int] = {}
    for line in bond_lines:
        left, right = (int(value) for value in line.split()[:2])
        counts[left] = counts.get(left, 0) + 1
        counts[right] = counts.get(right, 0) + 1
    return counts


def _angle_center_counts(angle_lines: list[str]) -> dict[int, int]:
    counts: dict[int, int] = {}
    for line in angle_lines:
        center = int(line.split()[1])
        counts[center] = counts.get(center, 0) + 1
    return counts


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
    dst = tmp_path / "system_itp" / "surface.itp"

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


def test_2_1_layers_scale_atom_count_and_z_levels(tmp_path):
    out = tmp_path / "surf_layers"

    dx = 0.47
    lx = 6.0
    ly = 6.0
    layers = 2
    dist_z = 0.4

    main([
        "--mode", "2-1",
        "--lx", str(lx),
        "--ly", str(ly),
        "--dx", str(dx),
        "--layers", str(layers),
        "--dist-z", str(dist_z),
        "--output", str(out),
    ])

    gro_file = tmp_path / "surf_layers.gro"
    lines = gro_file.read_text().splitlines()
    n_atoms = int(lines[1])

    scale = dx / 0.142
    a = 0.246 * scale
    lx_cell = 3.0 * a
    ly_cell = math.sqrt(3) * a
    nx = max(1, round(lx / lx_cell))
    ny = max(1, round(ly / ly_cell))
    expected = nx * ny * 6 * layers

    z_values = sorted({round(float(line[36:44]), 3) for line in lines[2:-1]})

    assert n_atoms == expected
    assert z_values == [3.0, 3.4]


def test_2_1_layers_default_to_sigma_based_spacing_for_regular_beads(tmp_path):
    out = tmp_path / "surf_sigma_regular"

    main([
        "--mode", "2-1",
        "--lx", "2",
        "--ly", "2",
        "--layers", "2",
        "--bead", "P4",
        "--output", str(out),
    ])

    lines = (tmp_path / "surf_sigma_regular.gro").read_text().splitlines()
    z_values = sorted({round(float(line[36:44]), 3) for line in lines[2:-1]})
    expected_spacing = round(0.47 * SIGMA_SPACING_FACTOR, 3)

    assert z_values == [3.0, round(3.0 + expected_spacing, 3)]


def test_4_1_layers_default_to_martini3_sigma_for_small_beads(tmp_path):
    out = tmp_path / "surf_sigma_small_m3"

    main([
        "--mode", "4-1",
        "--lx", "2",
        "--ly", "2",
        "--layers", "2",
        "--bead", "SC5",
        "--output", str(out),
    ])

    lines = (tmp_path / "surf_sigma_small_m3.gro").read_text().splitlines()
    z_values = sorted({round(float(line[36:44]), 3) for line in lines[2:-1]})
    expected_spacing = round(0.41 * SIGMA_SPACING_FACTOR, 3)

    assert z_values == [3.0, round(3.0 + expected_spacing, 3)]


def test_4_1_layers_can_use_martini2_sigma_for_dna_small_beads(tmp_path):
    out = tmp_path / "surf_sigma_small_m2"

    main([
        "--mode", "4-1",
        "--lx", "2",
        "--ly", "2",
        "--layers", "2",
        "--bead", "SC5",
        "--martini-version", "2",
        "--output", str(out),
    ])

    lines = (tmp_path / "surf_sigma_small_m2.gro").read_text().splitlines()
    z_values = sorted({round(float(line[36:44]), 3) for line in lines[2:-1]})
    expected_spacing = round(0.43 * SIGMA_SPACING_FACTOR, 3)

    assert z_values == [3.0, round(3.0 + expected_spacing, 3)]


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

    dst = tmp_path / "system_itp" / "surface.itp"
    assert dst.exists()


def test_graphene_periodic_mode(tmp_path):
    out = tmp_path / "graphene_surface"

    main([
        "--mode", "graphene",
        "--lx", "3",
        "--ly", "3",
        "--lz", "5",
        "--output", str(out),
    ])

    gro = tmp_path / "graphene_surface.gro"
    itp = tmp_path / "graphene_surface.itp"
    dst = tmp_path / "system_itp" / "surface.itp"

    assert gro.exists()
    assert itp.exists()
    assert dst.exists()
    assert "GRA" in gro.read_text()
    assert "[ bonds ]" in itp.read_text()


def test_local_multilayer_surface_cycles_beads_and_writes_full_itp(tmp_path):
    out = tmp_path / "layered_surface"

    main([
        "--mode", "4-1",
        "--lx", "3",
        "--ly", "3",
        "--layers", "2",
        "--bead", "P4", "C1",
        "--output", str(out),
    ])

    gro = tmp_path / "layered_surface.gro"
    itp = tmp_path / "layered_surface.itp"

    gro_text = gro.read_text()
    itp_text = itp.read_text()

    assert " P4" in gro_text
    assert " C1" in gro_text
    assert "P4" in itp_text
    assert "C1" in itp_text

    atom_count = int(gro_text.splitlines()[1].strip())
    atom_lines = _itp_section_lines(itp_text, "atoms")

    assert len(atom_lines) == atom_count


def test_local_uniform_surface_writes_expanded_itp_with_bonds_and_angles(tmp_path):
    out = tmp_path / "uniform_surface"

    main([
        "--mode", "4-1",
        "--lx", "3",
        "--ly", "3",
        "--bead", "C1",
        "--output", str(out),
    ])

    gro_text = (tmp_path / "uniform_surface.gro").read_text()
    itp_text = (tmp_path / "uniform_surface.itp").read_text()

    atom_count = int(gro_text.splitlines()[1].strip())
    atom_lines = _itp_section_lines(itp_text, "atoms")
    bond_lines = _itp_section_lines(itp_text, "bonds")
    angle_lines = _itp_section_lines(itp_text, "angles")

    assert len(atom_lines) == atom_count
    assert bond_lines
    assert angle_lines


def test_local_multilayer_surface_bonded_topology_stays_within_layers(tmp_path):
    out = tmp_path / "layered_surface_local_bonds"

    main([
        "--mode", "4-1",
        "--lx", "3",
        "--ly", "3",
        "--layers", "2",
        "--dist-z", "0.4",
        "--output", str(out),
    ])

    gro_text = (tmp_path / "layered_surface_local_bonds.gro").read_text()
    itp_text = (tmp_path / "layered_surface_local_bonds.itp").read_text()
    total_atoms = int(gro_text.splitlines()[1].strip())
    atoms_per_layer = total_atoms // 2

    for raw in _itp_section_lines(itp_text, "bonds"):
        left, right = (int(token) for token in raw.split()[:2])
        assert (left <= atoms_per_layer and right <= atoms_per_layer) or (
            left > atoms_per_layer and right > atoms_per_layer
        )

    for raw in _itp_section_lines(itp_text, "angles"):
        left, center, right = (int(token) for token in raw.split()[:3])
        assert (left <= atoms_per_layer and center <= atoms_per_layer and right <= atoms_per_layer) or (
            left > atoms_per_layer and center > atoms_per_layer and right > atoms_per_layer
        )


def test_local_periodic_xy_makes_edge_connectivity_uniform(tmp_path):
    out = tmp_path / "periodic_surface"

    main([
        "--mode", "4-1",
        "--lx", "3",
        "--ly", "3",
        "--output", str(out),
        "--periodic-xy",
    ])

    itp_text = (tmp_path / "periodic_surface.itp").read_text()
    bond_counts = _bond_degree_counts(_itp_section_lines(itp_text, "bonds"))
    angle_counts = _angle_center_counts(_itp_section_lines(itp_text, "angles"))

    assert bond_counts
    assert angle_counts
    assert set(bond_counts.values()) == {3}
    assert set(angle_counts.values()) == {3}


def test_2_1_layers_follow_graphite_like_abab_stacking(tmp_path):
    out = tmp_path / "surf_2_1_abab"

    dx = 0.47
    dist_z = 0.4
    main([
        "--mode", "2-1",
        "--lx", "1",
        "--ly", "1",
        "--dx", str(dx),
        "--layers", "3",
        "--dist-z", str(dist_z),
        "--output", str(out),
    ])

    lines = (tmp_path / "surf_2_1_abab.gro").read_text().splitlines()[2:-1]
    first_three = [lines[0], lines[6], lines[12]]
    coords = [(float(line[20:28]), float(line[28:36]), float(line[36:44])) for line in first_three]

    scale = dx / 0.142
    a = 0.246 * scale

    assert coords[0] == (0.0, 0.0, 3.0)
    assert coords[1] == (round(0.5 * a, 3), round((math.sqrt(3) / 6) * a, 3), 3.4)
    assert coords[2] == (0.0, 0.0, 3.8)


def test_4_1_layers_follow_graphite_like_abab_stacking(tmp_path):
    out = tmp_path / "surf_4_1_abab"

    dx = 0.47
    dist_z = 0.335
    main([
        "--mode", "4-1",
        "--lx", "1",
        "--ly", "1",
        "--dx", str(dx),
        "--layers", "3",
        "--dist-z", str(dist_z),
        "--output", str(out),
    ])

    lines = (tmp_path / "surf_4_1_abab.gro").read_text().splitlines()[2:-1]
    first_three = [lines[0], lines[4], lines[8]]
    coords = [(float(line[20:28]), float(line[28:36]), float(line[36:44])) for line in first_three]

    assert coords[0] == (0.0, 0.0, 3.0)
    assert coords[1] == (0.0, round(dx, 3), 3.335)
    assert coords[2] == (0.0, 0.0, 3.67)


def test_local_multilayer_surface_mixed_beads_use_average_sigma_between_layers(tmp_path):
    out = tmp_path / "surf_sigma_mixed"

    main([
        "--mode", "2-1",
        "--lx", "1",
        "--ly", "1",
        "--layers", "2",
        "--bead", "P4", "SC5",
        "--output", str(out),
    ])

    lines = (tmp_path / "surf_sigma_mixed.gro").read_text().splitlines()
    z_values = sorted({round(float(line[36:44]), 3) for line in lines[2:-1]})
    expected_spacing = round(((0.47 + 0.41) / 2.0) * SIGMA_SPACING_FACTOR, 3)

    assert z_values == [3.0, round(3.0 + expected_spacing, 3)]


def test_graphite_mode_copies_posre_support_files(tmp_path):
    out = tmp_path / "graphite_surface"

    main([
        "--mode", "graphite",
        "--lx", "3",
        "--ly", "3",
        "--graphite-layers", "3",
        "--graphite-spacing", "0.382",
        "--output", str(out),
    ])

    gro = tmp_path / "graphite_surface.gro"
    itp = tmp_path / "graphite_surface.itp"
    posre = tmp_path / "posre_GRA.itp"
    dst = tmp_path / "system_itp" / "surface.itp"
    dst_posre = tmp_path / "system_itp" / "posre_GRA.itp"

    assert gro.exists()
    assert itp.exists()
    assert posre.exists()
    assert dst.exists()
    assert dst_posre.exists()
    assert '#include "posre_GRA.itp"' in itp.read_text()


def test_cnt_m2_mode_generates_cnt_and_posres(tmp_path):
    out = tmp_path / "cnt_surface"

    main([
        "--mode", "cnt-m2",
        "--lx", "1",
        "--ly", "1",
        "--output", str(out),
    ])

    gro = tmp_path / "cnt_surface.gro"
    itp = tmp_path / "cnt_surface.itp"
    posres = tmp_path / "cnt_surface-posres.itp"
    dst = tmp_path / "system_itp" / "surface.itp"
    dst_posres = tmp_path / "system_itp" / "cnt_surface-posres.itp"

    assert gro.exists()
    assert itp.exists()
    assert posres.exists()
    assert dst.exists()
    assert dst_posres.exists()
    assert "CNP" in itp.read_text()
    assert '#include "cnt_surface-posres.itp"' in itp.read_text()


def test_cnt_m3_mode_uses_martini3_preset(tmp_path):
    out = tmp_path / "cnt_surface_m3"

    main([
        "--mode", "cnt-m3",
        "--lx", "1",
        "--ly", "1",
        "--output", str(out),
    ])

    gro = tmp_path / "cnt_surface_m3.gro"
    itp = tmp_path / "cnt_surface_m3.itp"

    assert gro.exists()
    assert itp.exists()
    assert "SC5" in itp.read_text()
    assert "SNda" not in itp.read_text()
    assert "0.410" in itp.read_text()
