import os
from pathlib import Path
import MDAnalysis as mda
import martinisurf.gromacs_inputs as gms


# ============================================================
# HELPERS
# ============================================================

def write_fake_gro(path: Path):
    """Create a minimal but valid GRO file for MDAnalysis."""
    gro = """Fake GRO
  2
    1ALA     CA    1   1.000   1.000   1.000
    2ALA     CA    2   2.000   2.000   2.000
   3.00000   3.00000   3.00000
"""
    path.write_text(gro)


def prepare_package_structure(tmp_path: Path):
    """
    Create the required package folders:
    surfmartini/
        ActiveITP/
        mdp_templates/
    """

    pkg = tmp_path / "surfmartini"
    pkg.mkdir()

    # Required to simulate the package
    (pkg / "__init__.py").write_text("")

    # Create ActiveITP with minimal dummy files
    A = pkg / "ActiveITP"
    A.mkdir()

    (A / "Active.itp").write_text("[ atoms ]\n1 ALA CA 1\n")
    (A / "martini_v3.0.0_Active.itp").write_text("dummy active")
    (A / "martini_v3.0.0_solvents_v1.itp").write_text("dummy solvent")
    (A / "martini_v3.0.0_ions_v1.itp").write_text("dummy ions")
    (A / "surface.itp").write_text("dummy surface")

    # Create mdp_templates
    mdp = pkg / "mdp_templates"
    mdp.mkdir()

    for name in ["nvt.mdp", "npt.mdp", "deposition.mdp", "production.mdp"]:
        (mdp / name).write_text("integrator = md\n")

    return pkg


# ============================================================
# TEST SUITE
# ============================================================

def test_gomartini_inside_2_system(tmp_path, monkeypatch):
    """
    Case 1 — running inside Simulation/2_system
    """
    pkg = prepare_package_structure(tmp_path)
    monkeypatch.setattr(gms, "surfmartini", __import__("surfmartini"))

    sim = tmp_path / "Simulation"
    sys2 = sim / "2_system"
    sys2.mkdir(parents=True)

    gro = sys2 / "Enzyme_Surface.gro"
    write_fake_gro(gro)

    # Switch working directory into Simulation/2_system
    monkeypatch.chdir(sys2)

    gms.main(["--anchor", "1", "2", "3"])

    assert (sim / "0_topology/system.top").exists()
    assert (sim / "0_topology/system_res.top").exists()
    assert (sim / "0_topology/index.ndx").exists()
    assert (sim / "1_mdp/nvt.mdp").exists()
    assert (sim / "2_system/system.gro").exists()


def test_gomartini_inside_Simulation(tmp_path, monkeypatch):
    """
    Case 2 — running inside Simulation/
    """
    pkg = prepare_package_structure(tmp_path)
    monkeypatch.setattr(gms, "surfmartini", __import__("surfmartini"))

    sim = tmp_path / "Simulation"
    sim.mkdir()

    sys2 = sim / "2_system"
    sys2.mkdir()

    gro = sys2 / "Enzyme_Surface.gro"
    write_fake_gro(gro)

    monkeypatch.chdir(sim)

    gms.main(["--anchor", "1", "2", "3"])

    assert (sim / "0_topology/system.top").exists()


def test_gomartini_with_outdir(tmp_path, monkeypatch):
    """
    Case 3 — running from outside using --outdir
    """
    pkg = prepare_package_structure(tmp_path)
    monkeypatch.setattr(gms, "surfmartini", __import__("surfmartini"))

    sim = tmp_path / "MySim"
    sys2 = sim / "2_system"
    sys2.mkdir(parents=True)

    gro = sys2 / "Enzyme_Surface.gro"
    write_fake_gro(gro)

    # Run from outside Simulation folder
    monkeypatch.chdir(tmp_path)

    gms.main(["--anchor", "1", "2", "3", "--outdir", str(sim)])

    assert (sim / "0_topology/system.top").exists()
    assert (sim / "0_topology/ActiveITP/Active_res.itp").exists()
    assert (sim / "1_mdp/deposition.mdp").exists()

