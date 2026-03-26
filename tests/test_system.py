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
    (itp / "martini_v3.0.0_solvents_v1.itp").write_text("dummy")
    (itp / "martini_v3.0.0_ions_v1.itp").write_text("dummy")

    mdp = sim / "mdp_templates"
    mdp.mkdir()
    for name in ["nvt.mdp", "npt.mdp", "deposition.mdp", "production.mdp"]:
        (mdp / name).write_text("integrator = md\n")

    return sim, sys2


def _normalize_mdp_nsteps(text: str) -> str:
    lines = []
    for line in text.splitlines():
        if line.strip().startswith("nsteps"):
            lines.append("nsteps                   = <ANY>")
        else:
            lines.append(line)
    return "\n".join(lines)


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


def test_water_template_is_copied_when_present(tmp_path, monkeypatch):
    sim, sys2 = prepare_simulation_structure(tmp_path)
    monkeypatch.chdir(sys2)

    # Add package-level water template expected by gromacs_inputs.
    pkg_templates = Path(gms.__file__).resolve().parent / "system_templates"
    pkg_templates.mkdir(exist_ok=True)
    water_tpl = pkg_templates / "water.gro"
    backup = water_tpl.read_text() if water_tpl.exists() else None
    water_tpl.write_text("WATER TEMPLATE TEST\n    0\n   1.00000   1.00000   1.00000\n")

    try:
        gms.main([
            "--moltype", "ENZ",
            "--anchor", "1", "2", "3",
        ])
        assert (sim / "2_system" / "water.gro").exists()
    finally:
        if backup is None:
            water_tpl.unlink(missing_ok=True)
        else:
            water_tpl.write_text(backup)


def test_polarizable_water_mode_uses_pw_templates_for_dna(tmp_path, monkeypatch):
    sim, sys2 = prepare_simulation_structure(tmp_path)
    monkeypatch.chdir(sys2)

    dna_gro = """DNA polarizable water test
  2
    1DG     BB1    1   1.000   1.000   1.000
    1DG     BB2    2   1.100   1.000   1.000
   3.00000   3.00000   3.00000
"""
    (sys2 / "immobilized_system.gro").write_text(dna_gro)
    (sim / "immobilized_system.gro").write_text(dna_gro)

    itp = sim / "system_itp"
    (itp / "Active.itp").unlink(missing_ok=True)
    (itp / "surface.itp").write_text("[ moleculetype ]\nSRF 1\n")
    (itp / "Nucleic_A+Nucleic_B.itp").write_text(
        "[ moleculetype ]\nNucleic_A+Nucleic_B 1\n\n"
        "[ atoms ]\n"
        "1 Q0 1 DG BB1 1 -1.0\n"
        "2 SN0 1 DG BB2 2 0.0\n"
    )
    (itp / "martini_v2.1-dna.itp").write_text("; standard dna\n")
    (itp / "martini_v2.1P-dna.itp").write_text("; polarizable dna\n")
    (itp / "martini_v2.0_ions.itp").write_text("; ions\n")

    gms.main([
        "--anchor", "1", "1",
        "--polarizable-water",
    ])

    assert (sim / "2_system" / "polarize-water.gro").exists()

    system_top = (sim / "0_topology" / "system.top").read_text()
    assert '#include "system_itp/martini_v2.1P-dna.itp"' in system_top
    assert '#include "system_itp/martini_v2.1-dna.itp"' not in system_top

    production = (sim / "1_mdp" / "production_dna.mdp").read_text()
    assert "coulombtype              = Cut-off" in production
    assert "coulomb-modifier         = Potential-shift" in production
    assert "vdwtype                  = Cut-off" in production
    assert "vdw-modifier             = Force-switch" in production
    assert "constraints              = none" in production
    assert "epsilon-r                = 2.5" in production
    assert "rvdw-switch              = 0.9" in production


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


def test_protein_without_linker_ignores_linker_specific_logic(tmp_path, monkeypatch):
    sim, _ = prepare_simulation_structure(tmp_path)
    monkeypatch.chdir(sim / "2_system")

    gms.main([
        "--moltype", "ENZ",
        "--anchor", "1", "1",
    ])

    index_text = (sim / "0_topology" / "index.ndx").read_text()
    system_top = (sim / "0_topology" / "system.top").read_text()

    assert "[ LINKER ]" not in index_text
    assert '#include "system_itp/linker.itp"' not in system_top


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


def test_molecules_block_adds_total_linkers_with_itp_moltype(tmp_path, monkeypatch):
    sim, _ = prepare_simulation_structure(tmp_path)
    monkeypatch.chdir(sim / "2_system")

    gro = """Linker count test
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
    (sim / "system_itp" / "peg_linker.itp").write_text(
        "[ moleculetype ]\nPEGX 1\n"
    )

    gms.main([
        "--moltype", "ENZ",
        "--anchor", "1", "1",
        "--anchor", "2", "2",
        "--use-linker",
        "--linker-resname", "LNK",
        "--linker-size", "2",
        "--linker-itp-name", "peg_linker.itp",
    ])

    system_top = (sim / "0_topology" / "system.top").read_text()
    system_res_top = (sim / "0_topology" / "system_res.top").read_text()

    assert "PEGX 2" in system_top
    assert "PEGX 2" in system_res_top


def test_linker_mode_restrains_linker_and_not_biomolecule(tmp_path, monkeypatch):
    sim, _ = prepare_simulation_structure(tmp_path)
    monkeypatch.chdir(sim / "2_system")

    gro = """Linker restraint routing test
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
    (sim / "system_itp" / "peg_linker.itp").write_text(
        "[ moleculetype ]\nPEGX 1\n\n"
        "[ atoms ]\n"
        "1 C1 1 LNK C1 1 0.0\n"
        "2 C1 1 LNK C2 2 0.0\n"
    )

    gms.main([
        "--moltype", "ENZ",
        "--anchor", "1", "1",
        "--anchor", "2", "2",
        "--use-linker",
        "--linker-resname", "LNK",
        "--linker-size", "2",
        "--linker-itp-name", "peg_linker.itp",
    ])

    system_res_top = (sim / "0_topology" / "system_res.top").read_text()
    linker_anchor_itp = (sim / "0_topology" / "system_itp" / "peg_linker_anchor.itp").read_text()

    assert '#include "system_itp/ENZ.itp"' in system_res_top
    assert '#include "system_itp/ENZ_anchor.itp"' not in system_res_top
    assert '#include "system_itp/peg_linker_anchor.itp"' in system_res_top
    assert "[ position_restraints ]" in linker_anchor_itp
    assert "2 1 1000 1000 0" in linker_anchor_itp


def test_molecules_block_falls_back_to_itp_stem_when_missing_moleculetype(tmp_path, monkeypatch):
    sim, _ = prepare_simulation_structure(tmp_path)
    monkeypatch.chdir(sim / "2_system")

    gro = """Linker count test fallback
  3
    1ALA     CA    1   1.000   1.000   1.200
    2LNK     C1    2   1.050   1.050   1.100
    2LNK     C2    3   1.050   1.050   0.300
   4.00000   4.00000   4.00000
"""
    (sim / "2_system" / "immobilized_system.gro").write_text(gro)
    (sim / "system_itp" / "peg_linker.itp").write_text("; dummy linker\n")

    gms.main([
        "--moltype", "ENZ",
        "--anchor", "1", "1",
        "--use-linker",
        "--linker-resname", "LNK",
        "--linker-size", "2",
        "--linker-itp-name", "peg_linker.itp",
    ])

    system_top = (sim / "0_topology" / "system.top").read_text()
    assert "peg_linker 1" in system_top


def test_substrate_itp_include_and_molecules_count(tmp_path, monkeypatch):
    sim, _ = prepare_simulation_structure(tmp_path)
    monkeypatch.chdir(sim / "2_system")

    (sim / "system_itp" / "substrate.itp").write_text(
        "[ moleculetype ]\nSUBX 1\n"
    )

    gms.main([
        "--moltype", "ENZ",
        "--anchor", "1", "1",
        "--substrate-itp-name", "substrate.itp",
        "--substrate-count", "3",
    ])

    system_top = (sim / "0_topology" / "system.top").read_text()
    system_res_top = (sim / "0_topology" / "system_res.top").read_text()

    assert '#include "system_itp/substrate.itp"' in system_top
    assert '#include "system_itp/substrate.itp"' in system_res_top
    assert "SUBX 3" in system_top
    assert "SUBX 3" in system_res_top


def test_substrate_molecules_fall_back_to_itp_stem(tmp_path, monkeypatch):
    sim, _ = prepare_simulation_structure(tmp_path)
    monkeypatch.chdir(sim / "2_system")

    (sim / "system_itp" / "cg_substrate.itp").write_text("; no moleculetype block\n")

    gms.main([
        "--moltype", "ENZ",
        "--anchor", "1", "1",
        "--substrate-itp-name", "cg_substrate.itp",
        "--substrate-count", "2",
    ])

    system_top = (sim / "0_topology" / "system.top").read_text()
    assert "cg_substrate 2" in system_top


def test_cofactor_itp_include_and_molecules_count(tmp_path, monkeypatch):
    sim, _ = prepare_simulation_structure(tmp_path)
    monkeypatch.chdir(sim / "2_system")

    (sim / "system_itp" / "NAD.itp").write_text(
        "[ moleculetype ]\nNAD 1\n"
    )

    gms.main([
        "--moltype", "ENZ",
        "--anchor", "1", "1",
        "--cofactor-itp-name", "NAD.itp",
        "--cofactor-count", "1",
    ])

    system_top = (sim / "0_topology" / "system.top").read_text()
    system_res_top = (sim / "0_topology" / "system_res.top").read_text()

    assert '#include "system_itp/NAD.itp"' in system_top
    assert '#include "system_itp/NAD.itp"' in system_res_top
    assert "NAD 1" in system_top
    assert "NAD 1" in system_res_top
    lines_top = [line.strip() for line in system_top.splitlines()]
    assert lines_top.index("ENZ 1") < lines_top.index("NAD 1")
    molecules_idx = lines_top.index("[ molecules ]")
    molecule_lines = [
        ln for ln in lines_top[molecules_idx + 1:]
        if ln and not ln.startswith(";")
    ]
    assert molecule_lines[0] == "ENZ 1"
    assert molecule_lines[1] == "NAD 1"


def test_cofactor_and_substrate_can_coexist(tmp_path, monkeypatch):
    sim, _ = prepare_simulation_structure(tmp_path)
    monkeypatch.chdir(sim / "2_system")

    (sim / "system_itp" / "NAD.itp").write_text(
        "[ moleculetype ]\nNAD 1\n"
    )
    (sim / "system_itp" / "substrate.itp").write_text(
        "[ moleculetype ]\nSUBX 1\n"
    )

    gms.main([
        "--moltype", "ENZ",
        "--anchor", "1", "1",
        "--cofactor-itp-name", "NAD.itp",
        "--cofactor-count", "1",
        "--substrate-itp-name", "substrate.itp",
        "--substrate-count", "3",
    ])

    system_top = (sim / "0_topology" / "system.top").read_text()
    assert "NAD 1" in system_top
    assert "SUBX 3" in system_top


def test_substrate_include_is_skipped_when_moltype_exists_in_forcefield(tmp_path, monkeypatch):
    sim, _ = prepare_simulation_structure(tmp_path)
    monkeypatch.chdir(sim / "2_system")

    # Simulate ETO already defined in the Martini solvent FF include.
    (sim / "system_itp" / "martini_v3.0.0_solvents_v1.itp").write_text(
        "[ moleculetype ]\nETO 1\n\n[ atoms ]\n1 SP1 1 ETO ETO 1 0.0\n"
    )
    (sim / "system_itp" / "ETO.itp").write_text(
        "[ moleculetype ]\nETO 1\n\n[ atoms ]\n1 SP1 1 ETO ETO 1 0.0\n"
    )

    gms.main([
        "--moltype", "ENZ",
        "--anchor", "1", "1",
        "--substrate-itp-name", "ETO.itp",
        "--substrate-count", "10",
    ])

    system_top = (sim / "0_topology" / "system.top").read_text()
    assert '#include "system_itp/ETO.itp"' not in system_top
    assert "ETO 10" in system_top


def test_substrate_conflicting_local_moltype_is_renamed_when_atomnames_differ(tmp_path, monkeypatch):
    sim, _ = prepare_simulation_structure(tmp_path)
    monkeypatch.chdir(sim / "2_system")

    (sim / "system_itp" / "martini_v3.0.0_solvents_v1.itp").write_text(
        "[ moleculetype ]\nETO 1\n\n[ atoms ]\n1 SP1 1 ETO ETO 1 0.0\n"
    )
    (sim / "system_itp" / "ETO.itp").write_text(
        "[ moleculetype ]\nETO 1\n\n[ atoms ]\n1 C1 1 ETO ETO1 1 0.0\n"
    )

    gms.main([
        "--moltype", "ENZ",
        "--anchor", "1", "1",
        "--substrate-itp-name", "ETO.itp",
        "--substrate-count", "10",
    ])

    system_top = (sim / "0_topology" / "system.top").read_text()
    local_itp = (sim / "0_topology" / "system_itp" / "ETO.itp").read_text()
    assert '#include "system_itp/ETO.itp"' in system_top
    assert "ETO_LOCAL 10" in system_top
    assert "ETO_LOCAL 1" in local_itp


def test_substrate_without_local_itp_uses_forcefield_moltype(tmp_path, monkeypatch):
    sim, _ = prepare_simulation_structure(tmp_path)
    monkeypatch.chdir(sim / "2_system")

    (sim / "system_itp" / "martini_v3.0.0_solvents_v1.itp").write_text(
        "[ moleculetype ]\nPPN 1\n"
    )

    gms.main([
        "--moltype", "ENZ",
        "--anchor", "1", "1",
        "--substrate-moltype", "PPN",
        "--substrate-count", "10",
    ])

    system_top = (sim / "0_topology" / "system.top").read_text()
    assert '#include "system_itp/PPN.itp"' not in system_top
    assert "PPN 10" in system_top


def test_restrained_topology_compatibility_alias_is_written(tmp_path, monkeypatch):
    sim, _ = prepare_simulation_structure(tmp_path)
    monkeypatch.chdir(sim / "2_system")

    gms.main([
        "--moltype", "ENZ",
        "--anchor", "1", "1",
    ])

    system_res_top = (sim / "0_topology" / "system_res.top").read_text()
    system_anchor_top = (sim / "0_topology" / "system_anchor.top").read_text()
    assert system_res_top == system_anchor_top


def test_legacy_gromacs_workflow_is_not_generated(tmp_path, monkeypatch):
    sim, _ = prepare_simulation_structure(tmp_path)
    monkeypatch.chdir(sim / "2_system")

    gms.main([
        "--moltype", "ENZ",
        "--anchor", "1", "1",
    ])

    assert not (sim / "1_mdp" / "gromacs_workflow").exists()


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
    assert "[ LNK ]" in index_text
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
    assert "pull-coord1-dim          = Y Y Y" in production
    assert "pull-coord2-dim          = N N Y" in production
    assert "pull-coord3-dim          = Y Y Y" in production
    assert "pull-coord4-dim          = N N Y" in production
    assert "[ Anchor_6 ]" in index_text


def test_linker_pull_uses_start_from_current_distance(tmp_path, monkeypatch):
    sim, _ = prepare_simulation_structure(tmp_path)
    monkeypatch.chdir(sim / "2_system")

    gro = """Linker start test
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
    assert "pull-coord1_start        = yes" in production
    assert "pull-coord2_start        = yes" in production
    assert "pull-coord1-dim          = Y Y Y" in production
    assert "pull-coord2-dim          = N N Y" in production
    assert "pull-coord1-init" not in production
    assert "pull-coord2-init" not in production


def test_linker_pull_writes_init_when_flags_are_provided(tmp_path, monkeypatch):
    sim, _ = prepare_simulation_structure(tmp_path)
    monkeypatch.chdir(sim / "2_system")

    gro = """Linker start test
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
        "--linker-pull-init-prot", "0.6",
        "--linker-pull-init-surf", "0.7",
    ])

    production = (sim / "1_mdp" / "production.mdp").read_text()
    deposition = (sim / "1_mdp" / "deposition.mdp").read_text()

    assert "pull-coord1-init         = 0.600" in production
    assert "pull-coord2-init         = 0.700" in production
    assert "pull-coord1-init         = 0.600" in deposition
    assert "pull-coord2-init         = 0.700" in deposition
    assert "pull-coord1_start        = yes" not in production
    assert "pull-coord2_start        = yes" not in production
    assert "pull-coord1_start        = yes" not in deposition
    assert "pull-coord2_start        = yes" not in deposition


def test_anchor_pull_can_include_init_distance():
    block = gms._build_pull_block(anchor_count=2, surface_group="SRF", init_nm=0.8)

    assert "pull-coord1_start        = yes" not in block
    assert "pull-coord2_start        = yes" not in block
    assert "pull-coord1-init         = 0.800" in block
    assert "pull-coord2-init         = 0.800" in block


def test_linker_pull_omits_surface_start_when_init_is_enabled():
    block_with_init = gms._build_linker_pull_block(
        linker_count=1,
        surface_group="SRF",
        init_prot_nm=0.6,
        init_surf_nm=0.7,
        include_init=True,
    )
    assert "pull-coord1_start        = yes" not in block_with_init
    assert "pull-coord2_start        = yes" not in block_with_init
    assert "pull-coord1-init         = 0.600" in block_with_init
    assert "pull-coord2-init         = 0.700" in block_with_init

    block_no_init = gms._build_linker_pull_block(
        linker_count=1,
        surface_group="SRF",
        init_prot_nm=0.6,
        init_surf_nm=0.7,
        include_init=False,
    )
    assert "pull-coord1_start        = yes" in block_no_init
    assert "pull-coord2_start        = yes" in block_no_init
    assert "pull-coord1-init" not in block_no_init
    assert "pull-coord2-init" not in block_no_init


def test_linker_pull_block_can_skip_head_coordinate_for_dna_bonded_mode():
    block = gms._build_linker_pull_block(
        linker_count=1,
        surface_group="SRF",
        init_prot_nm=0.6,
        init_surf_nm=0.7,
        include_init=True,
        include_head_pull=False,
    )

    assert "pull-ncoords             = 1" in block
    assert "pull-coord1-groups       = 3 4" in block
    assert "pull-coord1-dim          = N N Y" in block
    assert "pull-coord1-init         = 0.700" in block
    assert "pull-coord2-groups" not in block


def test_dna_linker_mode_uses_bonded_coupling_and_single_surface_pull(tmp_path, monkeypatch):
    sim, _ = prepare_simulation_structure(tmp_path)
    monkeypatch.chdir(sim / "2_system")

    gro = """DNA linker bonded test
  4
    1DG     BB1    1   1.000   1.000   1.000
    1DG     BB2    2   1.100   1.000   1.000
    2LNK     C1    3   1.020   1.000   0.900
    2LNK     C2    4   1.020   1.000   0.400
   3.00000   3.00000   3.00000
"""
    (sim / "2_system" / "immobilized_system.gro").write_text(gro)
    itp = sim / "system_itp"
    (itp / "Active.itp").unlink(missing_ok=True)
    (itp / "linker.itp").write_text(
        "[ moleculetype ]\nLNK 1\n\n"
        "[ atoms ]\n"
        "1 C1 1 LNK C1 1 0.0\n"
        "2 C1 1 LNK C2 2 0.0\n"
    )
    (itp / "Nucleic_A+Nucleic_B.itp").write_text(
        "[ moleculetype ]\nNucleic_A+Nucleic_B 1\n\n"
        "[ atoms ]\n"
        "1 Q0 1 DG BB1 1 -1.0\n"
        "2 SN0 1 DG BB2 2 0.0\n"
    )
    (itp / "martini_v2.1-dna.itp").write_text("; dummy\n")
    (itp / "martini_v2.0_ions.itp").write_text("; dummy\n")

    gms.main([
        "--anchor", "1", "1",
        "--use-linker",
        "--linker-resname", "LNK",
        "--linker-size", "2",
        "--linker-pull-init-prot", "0.47",
        "--linker-pull-init-surf", "0.60",
    ])

    production = (sim / "1_mdp" / "production_dna.mdp").read_text()
    assert "pull-ncoords             = 1" in production
    assert "pull-coord1-groups       = 3 4" in production
    assert "pull-coord1-dim          = N N Y" in production
    assert "pull-coord2-groups" not in production

    system_top = (sim / "0_topology" / "system.top").read_text()
    assert '#include "system_itp/Nucleic_A+Nucleic_B_linker.itp"' in system_top
    assert '#include "system_itp/linker.itp"' not in system_top
    assert "LNK 1" not in system_top

    merged_itp = sim / "0_topology" / "system_itp" / "Nucleic_A+Nucleic_B_linker.itp"
    text = merged_itp.read_text()
    assert "[ bonds ]" in text
    assert "[ angles ]" in text
    assert "[ position_restraints ]" in text
    assert "#ifdef POSRES_DNA" in text
    assert "4 1 1000 1000 0" in text
    assert "0.470 1250.0" in text
    assert "180.0 20.0" in text


def test_protein_deposition_and_production_templates_only_differ_in_nsteps():
    mdp_dir = Path(gms.__file__).resolve().parent / "mdp_templates"
    deposition = (mdp_dir / "deposition.mdp").read_text()
    production = (mdp_dir / "production.mdp").read_text()

    assert "compressibility          = 0 4e-5" in deposition
    assert "compressibility          = 0 4e-5" in production
    assert "define                   = -DPOSRES" not in deposition
    assert "define                   = -DPOSRES" in production

    normalized_deposition = _normalize_mdp_nsteps(deposition)
    normalized_production = _normalize_mdp_nsteps(production).replace("define                   = -DPOSRES\n\n", "", 1)
    assert normalized_deposition == normalized_production


def test_dna_deposition_and_production_templates_keep_expected_protocol_values():
    mdp_dir = Path(gms.__file__).resolve().parent / "mdp_templates"
    deposition = (mdp_dir / "deposition_dna.mdp").read_text()
    production = (mdp_dir / "production_dna.mdp").read_text()

    assert "dt                       = 0.01" in deposition
    assert "nsteps                   = 2000000" in deposition
    assert "dt                       = 0.01" in production
    assert "nsteps                   = 10000000" in production
    assert "coulombtype              = reaction-field" in deposition
    assert "coulombtype              = reaction-field" in production
    assert "vdw-modifier             = Potential-shift-verlet" in deposition
    assert "vdw-modifier             = Potential-shift-verlet" in production
    assert "tc-grps                  = System" in deposition
    assert "tc-grps                  = System" in production
    assert "ref-t                    = 300" in deposition
    assert "ref-t                    = 300" in production
    assert "Pcoupl                   = Parrinello-Rahman ;parrinello-rahman" in deposition
    assert "Pcoupltype               = semiisotropic" in deposition
    assert "tau-p                    = 12.0" in deposition
    assert "compressibility          = 0 3e-4" in deposition
    assert "ref-p                    = 1 1" in deposition
    assert "tau-p" not in production
    assert "compressibility" not in production
    assert "define                   = -DPOSRES" in deposition
    assert "define                   = -DPOSRES -DPOSRES_DNA" in production
    assert "freezegrps" not in deposition
    assert "freezegrps" not in production
    assert "freezedim" not in deposition
    assert "freezedim" not in production
    assert "cos-acceleration" not in deposition
    assert "cos-acceleration" not in production


def test_dna_nvt_template_uses_requested_thermostat_settings():
    mdp_dir = Path(gms.__file__).resolve().parent / "mdp_templates"
    nvt = (mdp_dir / "nvt_dna.mdp").read_text()

    assert "define                   = -DPOSRES" in nvt
    assert "POSRES_DNA" not in nvt
    assert "coulombtype              = reaction-field" in nvt
    assert "vdw-modifier             = Potential-shift-verlet" in nvt
    assert "tc-grps                  = System" in nvt
    assert "ref-t                    = 300" in nvt
    assert "freezegrps" not in nvt
    assert "freezedim" not in nvt


def test_dna_mdp_templates_omit_epsilon_rf():
    mdp_dir = Path(gms.__file__).resolve().parent / "mdp_templates"

    for name in (
        "minimization_dna.mdp",
        "nvt_dna.mdp",
        "npt_dna.mdp",
        "deposition_dna.mdp",
        "production_dna.mdp",
    ):
        text = (mdp_dir / name).read_text()
        assert "epsilon_rf" not in text


def test_dna_minimization_freezes_surface_and_npt_keeps_surface_posres(tmp_path):
    mdp_dir = Path(gms.__file__).resolve().parent / "mdp_templates"
    minimization = (mdp_dir / "minimization_dna.mdp").read_text()
    npt = (mdp_dir / "npt_dna.mdp").read_text()

    assert "define                   = -DPOSRES" not in minimization
    assert "define                   =" not in minimization
    assert "define                   = -DPOSRES" in npt
    assert "POSRES_DNA" not in minimization
    assert "POSRES_DNA" not in npt
    assert "freezegrps               = SRF" in minimization
    assert "freezedim                = Y Y Y" in minimization
    assert "freezegrps = SRF" in npt
    assert "freezedim = Y Y Y" in npt


def test_dna_surface_itp_gets_xyz_posres(tmp_path, monkeypatch):
    sim, sys2 = prepare_simulation_structure(tmp_path)
    monkeypatch.chdir(sys2)

    dna_gro = """DNA surface posres test
  5
    1DG     BB1    1   1.000   1.000   1.000
    1DG     BB2    2   1.100   1.000   1.000
    2SRF      C    3   1.500   1.500   0.500
    3SRF      C    4   2.000   1.500   0.500
    4SRF      C    5   2.500   1.500   0.500
   4.00000   4.00000   4.00000
"""
    (sys2 / "immobilized_system.gro").write_text(dna_gro)
    (sim / "immobilized_system.gro").write_text(dna_gro)

    itp = sim / "system_itp"
    (itp / "Active.itp").unlink(missing_ok=True)
    (itp / "surface.itp").write_text(
        "[ moleculetype ]\n"
        "SRF 1\n\n"
        "[ atoms ]\n"
        "1 C1 1 SRF C1 1 0.0\n"
        "2 C1 1 SRF C1 2 0.0\n"
        "3 C1 1 SRF C1 3 0.0\n\n"
        "[ bonds ]\n"
        "1 2 1 0.470 5000\n"
        "2 3 1 0.470 5000\n\n"
        "[ angles ]\n"
        "1 2 3 1 180.0 300\n"
    )
    (itp / "Nucleic_A+Nucleic_B.itp").write_text(
        "[ moleculetype ]\nNucleic_A+Nucleic_B 1\n\n"
        "[ atoms ]\n"
        "1 Q0 1 DG BB1 1 -1.0\n"
        "2 SN0 1 DG BB2 2 0.0\n"
    )
    (itp / "martini_v2.1-dna.itp").write_text("; standard dna\n")
    (itp / "martini_v2.1P-dna.itp").write_text("; polarizable dna\n")
    (itp / "martini_v2.0_ions.itp").write_text("; ions\n")

    gms.main([
        "--anchor", "1", "1",
    ])

    surface_itp = (sim / "0_topology" / "system_itp" / "surface.itp").read_text()
    dna_anchor_itp = (sim / "0_topology" / "system_itp" / "Nucleic_A+Nucleic_B_anchor.itp").read_text()
    deposition = (sim / "1_mdp" / "deposition_dna.mdp").read_text()
    production = (sim / "1_mdp" / "production_dna.mdp").read_text()
    system_top = (sim / "0_topology" / "system.top").read_text()

    assert "[ position_restraints ]" in surface_itp
    assert "#ifdef POSRES" in surface_itp
    assert "1 1 50000 50000 50000" in surface_itp
    assert "2 1 50000 50000 50000" in surface_itp
    assert "3 1 50000 50000 50000" in surface_itp
    assert "[ bonds ]" in surface_itp
    assert "[ angles ]" in surface_itp
    assert "[ position_restraints ]" in dna_anchor_itp
    assert "#ifdef POSRES_DNA" in dna_anchor_itp
    assert "SRF 1" in system_top
    assert "tc-grps                  = SRF DNA" in deposition
    assert "tau-t                    = 0.3 0.3" in deposition
    assert "ref-t                    = 300 300" in deposition
    assert "define                   = -DPOSRES" in deposition
    assert "POSRES_DNA" not in deposition
    assert "define                   = -DPOSRES -DPOSRES_DNA" in production
    assert "freezegrps" not in deposition


def test_dna_local_surface_itp_uses_uniform_50k_posres(tmp_path, monkeypatch):
    sim, sys2 = prepare_simulation_structure(tmp_path)
    monkeypatch.chdir(sys2)

    dna_gro = """DNA local surface posres test
  5
    1DG     BB1    1   1.000   1.000   1.000
    1DG     BB2    2   1.100   1.000   1.000
    2SRF      C    3   1.500   1.500   0.500
    3SRF      C    4   2.000   1.500   0.500
    4SRF      C    5   2.500   1.500   0.500
   4.00000   4.00000   4.00000
"""
    (sys2 / "immobilized_system.gro").write_text(dna_gro)
    (sim / "immobilized_system.gro").write_text(dna_gro)

    itp = sim / "system_itp"
    (itp / "Active.itp").unlink(missing_ok=True)
    (itp / "surface.itp").write_text(
        ";;;;;; Local bonded surface topology\n\n"
        "[ moleculetype ]\n"
        "SRF 1\n\n"
        "[ atoms ]\n"
        "1 C1 1 SRF C1 1 0.0\n"
        "2 C1 1 SRF C1 2 0.0\n"
        "3 C1 1 SRF C1 3 0.0\n\n"
        "[ bonds ]\n"
        "1 2 1 0.470 5000\n"
        "2 3 1 0.470 5000\n\n"
        "[ angles ]\n"
        "1 2 3 1 180.0 300\n"
    )
    (itp / "Nucleic_A+Nucleic_B.itp").write_text(
        "[ moleculetype ]\nNucleic_A+Nucleic_B 1\n\n"
        "[ atoms ]\n"
        "1 Q0 1 DG BB1 1 -1.0\n"
        "2 SN0 1 DG BB2 2 0.0\n"
    )
    (itp / "martini_v2.1-dna.itp").write_text("; standard dna\n")
    (itp / "martini_v2.1P-dna.itp").write_text("; polarizable dna\n")
    (itp / "martini_v2.0_ions.itp").write_text("; ions\n")

    gms.main([
        "--anchor", "1", "1",
    ])

    surface_itp = (sim / "0_topology" / "system_itp" / "surface.itp").read_text()

    assert "1 1 50000 50000 50000" in surface_itp
    assert "2 1 50000 50000 50000" in surface_itp
    assert "3 1 50000 50000 50000" in surface_itp


def test_refresh_dna_mdp_thermostat_groups_uses_topology_groups(tmp_path):
    mdp_dir = tmp_path / "1_mdp"
    itp_dir = tmp_path / "system_itp"
    top_path = tmp_path / "system_final.top"
    mdp_dir.mkdir()
    itp_dir.mkdir()

    for name in ("nvt_dna.mdp", "npt_dna.mdp", "deposition_dna.mdp", "production_dna.mdp"):
        (mdp_dir / name).write_text(
            "integrator = md\n"
            "rvdw                     = 1.1\n"
            "tcoupl                   = v-rescale\n"
            "tc-grps                  = System\n"
            "tau-t                    = 1.0\n"
            "ref-t                    = 300\n"
        )
    (mdp_dir / "minimization_dna.mdp").write_text(
        "integrator = steep\n"
        "define                   = -DPOSRES\n"
        "freezegrps               = SRF\n"
        "freezedim                = Y Y Y\n"
    )

    top_path.write_text(
        '#include "system_itp/Nucleic_A+Nucleic_B_linker.itp"\n'
        '#include "system_itp/surface.itp"\n\n'
        "[ molecules ]\n"
        "Nucleic_A+Nucleic_B 1\n"
        "SRF 1\n"
        "W 10\n"
        "WF 3\n"
        "NA 2\n"
        "NAD 1\n"
    )
    (itp_dir / "surface.itp").write_text(
        "[ moleculetype ]\nSRF 1\n\n[ atoms ]\n1 C1 1 SRF C1 1 0.0\n"
    )
    (itp_dir / "Nucleic_A+Nucleic_B.itp").write_text(
        "[ moleculetype ]\nNucleic_A+Nucleic_B 1\n\n"
        "[ atoms ]\n"
        "1 Q0 1 DG BB1 1 -1.0\n"
        "2 SN0 1 DG BB2 2 0.0\n"
    )
    (itp_dir / "Nucleic_A+Nucleic_B_linker.itp").write_text(
        "[ moleculetype ]\nNucleic_A+Nucleic_B 1\n\n"
        "[ atoms ]\n"
        "1 Q0 1 DG BB1 1 -1.0\n"
        "2 SN0 1 DG BB2 2 0.0\n"
        "3 C3 2 MOL1 C01 3 0.0\n"
    )
    (itp_dir / "NAD.itp").write_text(
        "[ moleculetype ]\nNAD 1\n\n"
        "[ atoms ]\n"
        "1 P5 1 NAD N1 1 0.0\n"
    )

    groups = gms.refresh_dna_mdp_thermostat_groups(
        mdp_dir=mdp_dir,
        top_path=top_path,
        itp_dir=itp_dir,
        surface_moltype="SRF",
        surface_resname="SRF",
    )

    production = (mdp_dir / "production_dna.mdp").read_text()
    minimization = (mdp_dir / "minimization_dna.mdp").read_text()

    assert groups == ["SRF", "DNA", "W", "WF", "IONS", "MOL1", "NAD"]
    assert "define                   = -DPOSRES" in (mdp_dir / "nvt_dna.mdp").read_text()
    assert "define                   = -DPOSRES -DPOSRES_DNA" in production
    assert "define                   = -DPOSRES" not in minimization
    assert "define                   =" not in minimization
    assert "freezegrps               = SRF" in minimization
    assert "freezedim                = Y Y Y" in minimization
    assert "tc-grps                  = SRF DNA W WF IONS MOL1 NAD" in production
    assert "tau-t                    = 0.3 0.3 0.3 0.3 0.3 0.3 0.3" in production
    assert "ref-t                    = 300 300 300 300 300 300 300" in production


def test_dna_generated_mdps_include_extra_topology_groups_and_index_aliases(tmp_path, monkeypatch):
    sim, sys2 = prepare_simulation_structure(tmp_path)
    monkeypatch.chdir(sys2)

    dna_gro = """DNA thermostat groups test
  7
    1DG     BB1    1   1.000   1.000   1.000
    1DG     BB2    2   1.100   1.000   1.000
    2SRF      C    3   1.500   1.500   0.500
    3SRF      C    4   2.000   1.500   0.500
    4SRF      C    5   2.500   1.500   0.500
    5MOL1   C01    6   1.800   2.200   1.400
    5MOL1   C02    7   1.900   2.200   1.100
   4.00000   4.00000   4.00000
"""
    (sys2 / "immobilized_system.gro").write_text(dna_gro)
    (sim / "immobilized_system.gro").write_text(dna_gro)

    itp = sim / "system_itp"
    (itp / "Active.itp").unlink(missing_ok=True)
    (itp / "surface.itp").write_text(
        "[ moleculetype ]\nSRF 1\n\n[ atoms ]\n1 C1 1 SRF C1 1 0.0\n"
    )
    (itp / "Nucleic_A+Nucleic_B.itp").write_text(
        "[ moleculetype ]\nNucleic_A+Nucleic_B 1\n\n"
        "[ atoms ]\n"
        "1 Q0 1 DG BB1 1 -1.0\n"
        "2 SN0 1 DG BB2 2 0.0\n"
        "3 C3 2 MOL1 C01 3 0.0\n"
        "4 C3 2 MOL1 C02 4 0.0\n"
    )
    (itp / "martini_v2.1-dna.itp").write_text("; standard dna\n")
    (itp / "martini_v2.1P-dna.itp").write_text("; polarizable dna\n")
    (itp / "martini_v2.0_ions.itp").write_text("; ions\n")

    gms.main([
        "--anchor", "1", "1",
    ])

    production = (sim / "1_mdp" / "production_dna.mdp").read_text()
    index_text = (sim / "0_topology" / "index.ndx").read_text()

    assert "tc-grps                  = SRF DNA MOL1" in production
    assert "tau-t                    = 0.3 0.3 0.3" in production
    assert "ref-t                    = 300 300 300" in production
    assert "[ MOL1 ]" in index_text


def test_write_custom_mdp_skips_pull_rewrite_when_disabled(tmp_path):
    src = tmp_path / "nvt_dna.mdp"
    dst = tmp_path / "nvt_dna_out.mdp"
    original = """integrator = md
pull                     = yes
pull-ngroups             = 2
pull-group1-name         = Anchor_1
pull-group2-name         = SRF
pull-ncoords             = 1
pull-coord1-geometry     = distance
pull-coord1-groups       = 1 2
pull-coord1-type         = umbrella
pull-coord1-k            = 1000.0
pull-coord1-rate         = 0
pull-coord1-dim          = N N Y
;pull-coord1_start        = yes
pull-coord1-init         = 0.8
"""
    src.write_text(original)

    gms.write_custom_mdp(
        src=src,
        dst=dst,
        anchor_count=1,
        is_dna=True,
        rewrite_pull=False,
    )

    assert dst.read_text() == original


def test_write_custom_mdp_can_write_init_when_enabled(tmp_path):
    src = tmp_path / "with_pull.mdp"
    dst = tmp_path / "with_pull_out.mdp"
    src.write_text(
        "integrator = md\n"
        "pull                     = yes\n"
        "pull-ngroups             = 2\n"
        "pull-group1-name         = Anchor_1\n"
        "pull-group2-name         = SRF\n"
        "pull-ncoords             = 1\n"
        "pull-coord1-geometry     = distance\n"
        "pull-coord1-groups       = 1 2\n"
        "pull-coord1-type         = umbrella\n"
        "pull-coord1-k            = 1000.0\n"
        "pull-coord1-rate         = 0\n"
        "pull-coord1-dim          = N N Y\n"
        "pull-coord1_start        = yes\n"
    )

    gms.write_custom_mdp(
        src=src,
        dst=dst,
        anchor_count=1,
        is_dna=True,
        anchor_pull_init_nm=0.8,
        rewrite_pull=True,
        include_pull_init=True,
    )

    out = dst.read_text()
    assert "pull-coord1_start        = yes" not in out
    assert "pull-coord1-init         = 0.800" in out


def test_materialize_posres_fc_in_itp_replaces_macro_with_numeric_value(tmp_path):
    itp = tmp_path / "Protein.itp"
    itp.write_text(
        "#ifdef POSRES\n"
        "#ifndef POSRES_FC\n"
        "#define POSRES_FC 750.0\n"
        "#endif\n"
        "[ position_restraints ]\n"
        "1 1 POSRES_FC POSRES_FC POSRES_FC\n"
        "#endif\n"
    )

    changed = gms._materialize_posres_fc_in_itp(itp)
    text = itp.read_text()

    assert changed is True
    assert "POSRES_FC POSRES_FC POSRES_FC" not in text
    assert "1 1 750.0 750.0 750.0" in text
    assert "#define POSRES_FC 750.0" in text


def test_rewrite_itp_with_posres_can_replace_legacy_block_and_macro(tmp_path):
    src = tmp_path / "Nucleic_A+Nucleic_B.itp"
    dst = tmp_path / "Nucleic_A+Nucleic_B_anchor.itp"
    src.write_text(
        "[ moleculetype ]\n"
        "Nucleic_A+Nucleic_B 1\n\n"
        "[ atoms ]\n"
        "1 Q0 1 DG BB1 1 -1.0\n"
        "2 SN0 1 DG BB2 2 0.0\n"
        "#ifdef POSRES\n"
        "#ifndef POSRES_FC\n"
        "#define POSRES_FC 750.0\n"
        "#endif\n"
        "[ position_restraints ]\n"
        "1 1 POSRES_FC POSRES_FC POSRES_FC\n"
        "#endif\n"
    )

    gms._rewrite_itp_with_posres(
        src_itp=src,
        dst_itp=dst,
        posres_atom_ids=[2],
        macro_name="POSRES_DNA",
    )

    text = dst.read_text()
    assert text.count("[ position_restraints ]") == 1
    assert "#ifdef POSRES_DNA" in text
    assert "#ifdef POSRES\n" not in text
    assert "2 1 1000 1000 0" in text


def test_dna_topologies_start_with_rubber_bands_define(tmp_path):
    topo_dir = tmp_path / "0_topology"
    itp_dir = topo_dir / "system_itp"
    topo_dir.mkdir(parents=True)
    itp_dir.mkdir(parents=True)

    # Minimal files required by write_top_files include discovery.
    for name in [
        "martini_v2.1-dna.itp",
        "martini_v2.1P-dna.itp",
        "martini_v2.0_ions.itp",
        "Nucleic_A+Nucleic_B.itp",
        "Nucleic_A+Nucleic_B_anchor.itp",
        "surface.itp",
    ]:
        (itp_dir / name).write_text("; dummy\n")

    gms.write_top_files(
        topo_dir=topo_dir,
        dst_itp_dir=itp_dir,
        moltype="Nucleic_A+Nucleic_B",
        mol_itp_name="Nucleic_A+Nucleic_B.itp",
        anchor_itp_name="Nucleic_A+Nucleic_B_anchor.itp",
        is_dna=True,
        use_linker=False,
    )

    system_top = (topo_dir / "system.top").read_text()
    system_res_top = (topo_dir / "system_res.top").read_text()

    assert system_top.startswith("#define RUBBER_BANDS")
    assert system_res_top.startswith("#define RUBBER_BANDS")


def test_protein_go_topologies_start_with_go_virt_define(tmp_path):
    topo_dir = tmp_path / "0_topology"
    itp_dir = topo_dir / "system_itp"
    topo_dir.mkdir(parents=True)
    itp_dir.mkdir(parents=True)

    for name in [
        "martini_v3.0.0.itp",
        "martini_v3.0.0_solvents_v1.itp",
        "martini_v3.0.0_ions_v1.itp",
        "Protein.itp",
        "Protein_anchor.itp",
        "surface.itp",
    ]:
        (itp_dir / name).write_text("; dummy\n")

    gms.write_top_files(
        topo_dir=topo_dir,
        dst_itp_dir=itp_dir,
        moltype="Protein",
        mol_itp_name="Protein.itp",
        anchor_itp_name="Protein_anchor.itp",
        is_dna=False,
        use_linker=False,
        go_model=True,
    )

    system_top = (topo_dir / "system.top").read_text()
    system_res_top = (topo_dir / "system_res.top").read_text()

    assert system_top.startswith("#define GO_VIRT")
    assert system_res_top.startswith("#define GO_VIRT")


def test_protein_topology_avoids_duplicate_defaults_includes(tmp_path):
    topo_dir = tmp_path / "0_topology"
    itp_dir = topo_dir / "system_itp"
    topo_dir.mkdir(parents=True)
    itp_dir.mkdir(parents=True)

    for name in [
        "martini_v3.0.0.itp",
        "martini_v3.0.0_solvents_v1.itp",
        "martini_v3.0.0_ions_v1.itp",
        "Protein.itp",
        "Protein_anchor.itp",
        "surface.itp",
    ]:
        (itp_dir / name).write_text("; dummy\n")

    gms.write_top_files(
        topo_dir=topo_dir,
        dst_itp_dir=itp_dir,
        moltype="Protein",
        mol_itp_name="Protein.itp",
        anchor_itp_name="Protein_anchor.itp",
        is_dna=False,
        use_linker=False,
    )

    system_top = (topo_dir / "system.top").read_text()
    assert '#include "system_itp/martini_v3.0.0.itp"' in system_top
    assert '#include "system_itp/martini_v3.0.0_Active.itp"' not in system_top


def test_surface_molecule_count_tracks_surface_atoms_for_single_atom_surface_itp(tmp_path, monkeypatch):
    sim, _ = prepare_simulation_structure(tmp_path)
    monkeypatch.chdir(sim / "2_system")

    gro = """Surface count test
  4
    1ALA     CA    1   1.000   1.000   1.000
    2SRF      C    2   1.500   1.500   0.500
    3SRF      C    3   2.000   1.500   0.500
    4SRF      C    4   2.500   1.500   0.500
   4.00000   4.00000   4.00000
"""
    (sim / "2_system" / "immobilized_system.gro").write_text(gro)
    (sim / "system_itp" / "surface.itp").write_text(
        "[ moleculetype ]\nSRF 1\n\n[ atoms ]\n1 C1 1 SRF C 1 0.0\n"
    )
    (sim / "system_itp" / "Protein_0.itp").write_text(
        "[ moleculetype ]\nPROT_REAL 1\n\n[ atoms ]\n1 C1 1 PRO A1 1 0.0\n"
    )

    gms.main([
        "--moltype", "Protein",
        "--anchor", "1", "1",
    ])

    system_top = (sim / "0_topology" / "system.top").read_text()
    assert "SRF 3" in system_top


def test_surface_molecule_count_ignores_legacy_tilde_lines_in_surface_itp(tmp_path, monkeypatch):
    sim, _ = prepare_simulation_structure(tmp_path)
    monkeypatch.chdir(sim / "2_system")

    gro = """Surface count test
  4
    1ALA     CA    1   1.000   1.000   1.000
    2SRF      C    2   1.500   1.500   0.500
    3SRF      C    3   2.000   1.500   0.500
    4SRF      C    4   2.500   1.500   0.500
   4.00000   4.00000   4.00000
"""
    (sim / "2_system" / "immobilized_system.gro").write_text(gro)
    (sim / "system_itp" / "surface.itp").write_text(
        "[ moleculetype ]\nSRF 1\n\n[ atoms ]\n1 C1 1 SRF C 1 0.0\n\n~\n"
    )
    (sim / "system_itp" / "Protein_0.itp").write_text(
        "[ moleculetype ]\nPROT_REAL 1\n\n[ atoms ]\n1 C1 1 PRO A1 1 0.0\n"
    )

    gms.main([
        "--moltype", "Protein",
        "--anchor", "1", "1",
    ])

    system_top = (sim / "0_topology" / "system.top").read_text()
    assert "SRF 3" in system_top


def test_non_srf_surface_gets_srf_alias_in_index(tmp_path, monkeypatch):
    sim, _ = prepare_simulation_structure(tmp_path)
    monkeypatch.chdir(sim / "2_system")

    gro = """CNT alias test
  4
    1ALA     CA    1   1.000   1.000   1.000
    2CNT      C    2   1.500   1.500   0.500
    3CNT      C    3   2.000   1.500   0.500
    4CNT      C    4   2.500   1.500   0.500
   4.00000   4.00000   4.00000
"""
    (sim / "2_system" / "immobilized_system.gro").write_text(gro)
    (sim / "system_itp" / "surface.itp").write_text(
        "[ moleculetype ]\n"
        "cnt-24-9-f11-SC5 1\n\n"
        "[ atoms ]\n"
        "1 SC5 1 CNT C 1 0.0\n"
    )
    (sim / "system_itp" / "Protein_0.itp").write_text(
        "[ moleculetype ]\nPROT_REAL 1\n\n[ atoms ]\n1 C1 1 PRO A1 1 0.0\n"
    )

    gms.main([
        "--moltype", "Protein",
        "--anchor", "1", "1",
    ])

    index_text = (sim / "0_topology" / "index.ndx").read_text()
    assert "[ cnt-24-9-f11-SC5 ]" in index_text
    assert "[ SRF ]" in index_text


def test_missing_named_protein_itp_uses_real_candidate_but_keeps_requested_moltype(tmp_path, monkeypatch):
    sim, _ = prepare_simulation_structure(tmp_path)
    monkeypatch.chdir(sim / "2_system")

    # Force missing Protein.itp path and provide a realistic martinize output instead.
    (sim / "system_itp" / "Active.itp").unlink(missing_ok=True)
    (sim / "system_itp" / "Protein_0.itp").write_text(
        "[ moleculetype ]\n"
        "Protein_0 1\n\n"
        "[ atoms ]\n"
        "1 C1 1 PRO A1 1 0.0\n"
        "2 C1 1 PRO A2 2 0.0\n"
        "3 C1 1 PRO A3 3 0.0\n"
    )

    gms.main([
        "--moltype", "Protein",
        "--anchor", "1", "1",
    ])

    protein_itp = sim / "0_topology" / "system_itp" / "Protein.itp"
    assert protein_itp.exists()
    text = protein_itp.read_text()
    assert "Protein 1" in text
    assert "Protein_0 1" not in text
