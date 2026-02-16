from pathlib import Path
import martinisurf.gromacs_inputs as gms


# ============================================================
# HELPERS
# ============================================================

def write_fake_gro(path: Path):
    """Create a minimal but valid GRO file."""
    gro = """Fake GRO
  2
    1ALA     CA    1   1.000   1.000   1.000
    2ALA     CA    2   2.000   2.000   2.000
   3.00000   3.00000   3.00000
"""
    path.write_text(gro)


def prepare_simulation_structure(base: Path):
    sim = base / "Simulation"
    sys2 = sim / "2_system"
    sys2.mkdir(parents=True)

    # CLI expects THIS exact name
    write_fake_gro(sys2 / "immobilized_system.gro")
    write_fake_gro(sim / "immobilized_system.gro")

    # Required folders
    (sim / "0_topology").mkdir()
    (sim / "1_mdp").mkdir()
    (sim / "system_itp").mkdir()

    itp = sim / "system_itp"
    (itp / "surface.itp").write_text("dummy surface")
    (itp / "Active.itp").write_text("dummy active")
    (itp / "martini_v3.0.0_Active.itp").write_text("dummy")
    (itp / "martini_v3.0.0_solvents_v1.itp").write_text("dummy")
    (itp / "martini_v3.0.0_ions_v1.itp").write_text("dummy")

    mdp = sim / "mdp_templates"
    mdp.mkdir()
    for name in ["nvt.mdp", "npt.mdp", "deposition.mdp", "production.mdp"]:
        (mdp / name).write_text("integrator = md\n")

    return sim, sys2


# ============================================================
# TEST SUITE
# ============================================================

def test_gomartini_inside_2_system(tmp_path, monkeypatch):
    sim, sys2 = prepare_simulation_structure(tmp_path)

    monkeypatch.chdir(sys2)

    gms.main([
        "--moltype", "ENZ",
        "--anchor", "1", "2", "3"
    ])

    assert (sim / "0_topology/system.top").exists()
    assert (sim / "0_topology/system_res.top").exists()
    assert (sim / "0_topology/index.ndx").exists()
    assert (sim / "1_mdp/nvt.mdp").exists()

    # IMPORTANT: output is written relative to cwd
    assert (sys2 / "system.gro").exists()


def test_gomartini_inside_Simulation(tmp_path, monkeypatch):
    sim, sys2 = prepare_simulation_structure(tmp_path)

    monkeypatch.chdir(sim)

    gms.main([
        "--moltype", "ENZ",
        "--anchor", "1", "2", "3"
    ])

    assert (sim / "0_topology/system.top").exists()
    assert (sim / "2_system/system.gro").exists()


def test_gomartini_with_outdir(tmp_path, monkeypatch):
    sim, sys2 = prepare_simulation_structure(tmp_path)

    monkeypatch.chdir(tmp_path)

    gms.main([
        "--moltype", "ENZ",
        "--anchor", "1", "2", "3",
        "--outdir", str(sim)
    ])

    assert (sim / "0_topology/system.top").exists()
    assert (sim / "1_mdp/deposition.mdp").exists()
    assert (sim / "2_system/system.gro").exists()


def test_dynamic_pull_single_anchor(tmp_path, monkeypatch):
    sim, _ = prepare_simulation_structure(tmp_path)
    monkeypatch.chdir(sim / "2_system")

    gms.main([
        "--moltype", "ENZ",
        "--anchor", "1", "1",
    ])

    production = (sim / "1_mdp" / "production.mdp").read_text()
    assert "pull-ncoords             = 1" in production
    assert "pull-coord1-groups       = 1 2" in production
    assert "pull-coord2-groups" not in production


def test_linker_itp_include_when_enabled(tmp_path, monkeypatch):
    sim, _ = prepare_simulation_structure(tmp_path)
    monkeypatch.chdir(sim / "2_system")

    topo_itp = sim / "system_itp"
    (topo_itp / "linker.itp").write_text("; dummy linker\n")

    gms.main([
        "--moltype", "ENZ",
        "--anchor", "1", "1",
        "--anchor", "2", "2",
        "--use-linker",
    ])

    system_top = (sim / "0_topology" / "system.top").read_text()
    system_res_top = (sim / "0_topology" / "system_res.top").read_text()

    assert '#include "system_itp/linker.itp"' in system_top
    assert '#include "system_itp/linker.itp"' in system_res_top


def test_linker_itp_include_uses_cli_name(tmp_path, monkeypatch):
    sim, _ = prepare_simulation_structure(tmp_path)
    monkeypatch.chdir(sim / "2_system")

    topo_itp = sim / "system_itp"
    (topo_itp / "peg_linker.itp").write_text("; dummy linker\n")

    gms.main([
        "--moltype", "ENZ",
        "--anchor", "1", "1",
        "--use-linker",
        "--linker-itp-name", "peg_linker.itp",
    ])

    system_top = (sim / "0_topology" / "system.top").read_text()
    system_res_top = (sim / "0_topology" / "system_res.top").read_text()

    assert '#include "system_itp/peg_linker.itp"' in system_top
    assert '#include "system_itp/peg_linker.itp"' in system_res_top


def test_linker_pull_generates_two_coordinates(tmp_path, monkeypatch):
    sim, _ = prepare_simulation_structure(tmp_path)
    monkeypatch.chdir(sim / "2_system")

    # One selected biomolecule residue (ALA) + one linker residue (LNK)
    gro = """Linker pull test
  2
    1ALA     CA    1   1.000   1.000   1.000
    2LNK     C1    2   1.500   1.500   0.800
   3.00000   3.00000   3.00000
"""
    (sim / "2_system" / "immobilized_system.gro").write_text(gro)
    (sim / "system_itp" / "linker.itp").write_text("; dummy linker\n")

    gms.main([
        "--moltype", "ENZ",
        "--anchor", "1", "1",
        "--use-linker",
        "--linker-resname", "LNK",
        "--linker-size", "1",
    ])

    production = (sim / "1_mdp" / "production.mdp").read_text()
    index_text = (sim / "0_topology" / "index.ndx").read_text()

    assert "pull-ncoords             = 2" in production
    assert "pull-coord1-groups       = 1 2" in production
    assert "pull-coord2-groups       = 3 4" in production
    assert "[ Anchor_2 ]" in index_text
    assert "[ Anchor_3 ]" in index_text
    # Anchor_1 must remain the selected biomolecule atoms only (not linker atom 2).
    anchor1_block = index_text.split("[ Anchor_1 ]", 1)[1].split("[ Anchor_2 ]", 1)[0]
    assert "1" in anchor1_block
    assert "2" not in anchor1_block


def test_multi_linker_generates_3n_anchors_and_2n_pulls(tmp_path, monkeypatch):
    sim, _ = prepare_simulation_structure(tmp_path)
    monkeypatch.chdir(sim / "2_system")

    # Two selected residues (ALA 1, ALA 2) and two linker instances (2 atoms each).
    gro = """Multi linker pull test
  6
    1ALA     CA    1   1.000   1.000   1.200
    2ALA     CA    2   2.000   2.000   1.200
    3LNK     C1    3   1.050   1.050   1.100
    3LNK     C2    4   1.050   1.050   0.300
    4LNK     C1    5   2.050   2.050   1.100
    4LNK     C2    6   2.050   2.050   0.300
   4.00000   4.00000   4.00000
"""
    (sim / "2_system" / "immobilized_system.gro").write_text(gro)
    (sim / "system_itp" / "linker.itp").write_text("; dummy linker\n")

    gms.main([
        "--moltype", "ENZ",
        "--anchor", "1", "1",
        "--anchor", "2", "2",
        "--use-linker",
        "--linker-resname", "LNK",
        "--linker-size", "2",
    ])

    production = (sim / "1_mdp" / "production.mdp").read_text()
    index_text = (sim / "0_topology" / "index.ndx").read_text()

    assert "pull-ngroups             = 7" in production  # 3N + 1 with N=2
    assert "pull-ncoords             = 4" in production  # 2N with N=2
    assert "pull-coord1-groups       = 1 2" in production
    assert "pull-coord2-groups       = 3 7" in production
    assert "pull-coord3-groups       = 4 5" in production
    assert "pull-coord4-groups       = 6 7" in production
    assert "[ Anchor_6 ]" in index_text


def test_linker_pull_init_values_are_configurable(tmp_path, monkeypatch):
    sim, _ = prepare_simulation_structure(tmp_path)
    monkeypatch.chdir(sim / "2_system")

    gro = """Linker init test
  2
    1ALA     CA    1   1.000   1.000   1.000
    2LNK     C1    2   1.500   1.500   0.800
   3.00000   3.00000   3.00000
"""
    (sim / "2_system" / "immobilized_system.gro").write_text(gro)
    (sim / "system_itp" / "linker.itp").write_text("; dummy linker\n")

    gms.main([
        "--moltype", "ENZ",
        "--anchor", "1", "1",
        "--use-linker",
        "--linker-resname", "LNK",
        "--linker-size", "1",
        "--linker-pull-init-prot", "0.430",
        "--linker-pull-init-surf", "0.365",
    ])

    production = (sim / "1_mdp" / "production.mdp").read_text()
    assert "pull-coord1-init         = 0.430" in production
    assert "pull-coord2-init         = 0.365" in production
