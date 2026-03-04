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
        "[ moleculetype ]\nETO 1\n"
    )
    (sim / "system_itp" / "ETO.itp").write_text(
        "[ moleculetype ]\nETO 1\n"
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


def test_linker_pull_ignores_legacy_init_flags(tmp_path, monkeypatch):
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
    assert "pull-coord1_start        = yes" in production
    assert "pull-coord2_start        = yes" in production
    assert "pull-coord1-init" not in production
    assert "pull-coord2-init" not in production


def test_anchor_pull_always_uses_start_from_current_distance():
    block = gms._build_pull_block(anchor_count=2, surface_group="SRF", init_nm=0.8)

    assert "pull-coord1_start        = yes" in block
    assert "pull-coord2_start        = yes" in block
    assert "pull-coord1-init" not in block
    assert "pull-coord2-init" not in block


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


def test_dna_deposition_and_production_templates_only_differ_in_nsteps():
    mdp_dir = Path(gms.__file__).resolve().parent / "mdp_templates"
    deposition = (mdp_dir / "deposition_dna.mdp").read_text()
    production = (mdp_dir / "production_dna.mdp").read_text()

    assert "compressibility          = 0 4e-5" in deposition
    assert "compressibility          = 0 4e-5" in production
    assert "define                   = -DPOSRES" not in deposition
    assert "define                   = -DPOSRES" in production

    normalized_deposition = _normalize_mdp_nsteps(deposition)
    normalized_production = _normalize_mdp_nsteps(production).replace("define                   = -DPOSRES\n\n", "", 1)
    assert normalized_deposition == normalized_production


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
