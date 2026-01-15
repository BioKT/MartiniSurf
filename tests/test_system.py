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
    """
    Create a minimal Simulation layout expected by gromacs_inputs.
    """
    sim = base / "Simulation"
    sys2 = sim / "2_system"
    sys2.mkdir(parents=True)

    # Fake enzyme–surface structure
    gro = sys2 / "Enzyme_Surface.gro"
    write_fake_gro(gro)

    # Required folders
    (sim / "0_topology").mkdir()
    (sim / "1_mdp").mkdir()
    (sim / "system_itp").mkdir()

    # Minimal ITPs
    itp = sim / "system_itp"
    (itp / "surface.itp").write_text("dummy surface")
    (itp / "Active.itp").write_text("dummy active")
    (itp / "martini_v3.0.0_Active.itp").write_text("dummy")
    (itp / "martini_v3.0.0_solvents_v1.itp").write_text("dummy")
    (itp / "martini_v3.0.0_ions_v1.itp").write_text("dummy")

    # MDP templates
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
    assert (sim / "2_system/system.gro").exists()


def test_gomartini_inside_Simulation(tmp_path, monkeypatch):
    sim, sys2 = prepare_simulation_structure(tmp_path)

    monkeypatch.chdir(sim)

    gms.main([
        "--moltype", "ENZ",
        "--anchor", "1", "2", "3"
    ])

    assert (sim / "0_topology/system.top").exists()


def test_gomartini_with_outdir(tmp_path, monkeypatch):
    sim, sys2 = prepare_simulation_structure(tmp_path)

    monkeypatch.chdir(tmp_path)

    gms.main([
        "--moltype", "ENZ",
        "--anchor", "1", "2", "3",
        "--outdir", str(sim)
    ])

    assert (sim / "0_topology/system.top").exists()
    assert (sim / "system_itp/Active_res.itp").exists()
    assert (sim / "1_mdp/deposition.mdp").exists()

