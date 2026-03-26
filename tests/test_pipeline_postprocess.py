import subprocess
from argparse import Namespace
from pathlib import Path

import pytest

from martinisurf import pipeline


def test_parser_accepts_solvate_flags():
    parser = pipeline.build_parser()
    args = parser.parse_args([
        "--pdb", "1RJW",
        "--anchor", "1", "1",
        "--lx", "10",
        "--ly", "10",
        "--solvate",
        "--salt-conc", "0.20",
        "--solvate-radius", "0.21",
    ])
    pipeline._validate_args(parser, args)
    assert args.solvate is True
    assert args.ionize is False
    assert args.salt_conc == 0.20
    assert args.solvate_radius == 0.21
    assert args.solvate_surface_clearance == 0.4


def test_parser_uses_dna_default_solvate_surface_clearance():
    parser = pipeline.build_parser()
    args = parser.parse_args([
        "--dna",
        "--pdb", "4C64",
        "--anchor", "1", "1",
        "--lx", "10",
        "--ly", "10",
        "--solvate",
    ])
    pipeline._validate_args(parser, args)
    assert args.solvate_surface_clearance == 0.4


def test_parser_keeps_explicit_dna_solvate_surface_clearance():
    parser = pipeline.build_parser()
    args = parser.parse_args([
        "--dna",
        "--pdb", "4C64",
        "--anchor", "1", "1",
        "--lx", "10",
        "--ly", "10",
        "--solvate",
        "--solvate-surface-clearance", "0.2",
    ])
    pipeline._validate_args(parser, args)
    assert args.solvate_surface_clearance == 0.2


def test_parser_accepts_complex_config_without_pdb(tmp_path):
    parser = pipeline.build_parser()
    surface = tmp_path / "surface.gro"
    surface.write_text("surface\n0\n1 1 1\n")
    cfg = tmp_path / "complex_config.yaml"
    cfg.write_text("mode: pre_cg_complex\n")
    args = parser.parse_args([
        "--complex-config", str(cfg),
        "--surface", str(surface),
    ])
    pipeline._validate_args(parser, args)


def test_parser_accepts_complex_config_with_linker_mode(tmp_path):
    parser = pipeline.build_parser()
    surface = tmp_path / "surface.gro"
    surface.write_text("surface\n0\n1 1 1\n")
    cfg = tmp_path / "complex_config.yaml"
    cfg.write_text("mode: pre_cg_complex\n")
    args = parser.parse_args([
        "--complex-config", str(cfg),
        "--surface", str(surface),
        "--linker", "input/linker.gro",
        "--linker-group", "1", "8", "10", "11",
    ])
    pipeline._validate_args(parser, args)


def test_preconfig_balance_low_z_only_applies_in_pure_preconfig_orientation():
    complex_cfg = {"balance_low_z": True}

    args = Namespace(linker=None, anchor=None, ads_mode=False)
    assert pipeline._use_preconfig_balance_low_z(args, complex_cfg) is True

    args = Namespace(linker="input/linker.gro", anchor=None, ads_mode=False)
    assert pipeline._use_preconfig_balance_low_z(args, complex_cfg) is False

    args = Namespace(linker=None, anchor=[["1", "8", "10"]], ads_mode=False)
    assert pipeline._use_preconfig_balance_low_z(args, complex_cfg) is False

    args = Namespace(linker=None, anchor=None, ads_mode=True)
    assert pipeline._use_preconfig_balance_low_z(args, complex_cfg) is False


def test_anchor_landmark_mode_uses_config_only_for_preconfig():
    args = Namespace(linker=None, anchor=[["A", "8", "10", "11"]], ads_mode=False)
    assert pipeline._anchor_landmark_mode_for_pipeline(args, None) is None

    args = Namespace(linker=None, anchor=None, ads_mode=False)
    assert pipeline._anchor_landmark_mode_for_pipeline(args, {"anchor_landmark_mode": "group"}) == "group"
    assert pipeline._anchor_landmark_mode_for_pipeline(args, {"anchor_landmark_mode": "residue"}) == "residue"


def test_parser_accepts_chain_based_anchor_and_linker_groups():
    parser = pipeline.build_parser()
    args = parser.parse_args([
        "--pdb", "1RJW",
        "--anchor", "B", "8", "10", "11",
        "--linker", "input/linker.gro",
        "--linker-group", "D", "8", "10", "11",
        "--lx", "10",
        "--ly", "10",
    ])

    assert args.anchor == [["B", "8", "10", "11"]]
    assert args.linker_group == [["D", "8", "10", "11"]]


def test_pipeline_reads_substrate_moltype_from_gro_when_itp_is_missing(tmp_path):
    gro = tmp_path / "PPN.gro"
    gro.write_text(
        "PPN\n"
        "    1\n"
        "    1PPN    C1    1   0.100   0.100   0.100\n"
        "   1.00000   1.00000   1.00000\n"
    )

    assert pipeline._read_gro_first_resname(str(gro)) == "PPN"


def test_normalize_uniform_atom_names_uses_top_includes_not_unincluded_itps(tmp_path):
    top_dir = tmp_path / "0_topology"
    itp_dir = top_dir / "system_itp"
    itp_dir.mkdir(parents=True)
    top = top_dir / "system_final.top"
    gro = tmp_path / "system_final.gro"

    (itp_dir / "martini_v3.0.0_solvents_v1.itp").write_text(
        "[ moleculetype ]\nETO 1\n\n[ atoms ]\n1 SP1 1 ETO ETO 1 0.0\n"
    )
    (itp_dir / "ETO.itp").write_text(
        "[ moleculetype ]\nETO 1\n\n[ atoms ]\n1 C1 1 ETO ETO1 1 0.0\n"
    )
    top.write_text(
        '#include "system_itp/martini_v3.0.0_solvents_v1.itp"\n\n'
        "[ system ]\nX\n\n"
        "[ molecules ]\n"
        "ETO 1\n"
    )
    gro.write_text(
        "ETO test\n"
        "    1\n"
        "    1ETO   ETO1    1   0.100   0.100   0.100\n"
        "   1.00000   1.00000   1.00000\n"
    )

    changed = pipeline._normalize_uniform_atom_names_from_itp(top_dir=top_dir, top_path=top, gro_path=gro)

    assert changed is True
    assert "ETO1" not in gro.read_text()
    assert "ETO" in gro.read_text()


def test_normalize_surface_itp_atomnames_matches_uniform_surface_gro_name(tmp_path):
    surface_gro = tmp_path / "surface.gro"
    surface_itp = tmp_path / "surface.itp"
    surface_gro.write_text(
        "Surface\n"
        "    2\n"
        "    1GRA      C    1   0.100   0.100   0.100\n"
        "    2GRA      C    2   0.200   0.200   0.100\n"
        "   1.00000   1.00000   1.00000\n"
    )
    surface_itp.write_text(
        "[ moleculetype ]\n"
        "GRA 1\n\n"
        "[ atoms ]\n"
        "1 C1 1 GRA C1 1 0.0\n"
    )

    changed = pipeline._normalize_surface_itp_atomnames(surface_gro=surface_gro, surface_itp=surface_itp)

    assert changed is True
    assert "1 C1 1 GRA C 1 0.0" in surface_itp.read_text()


def test_normalize_surface_itp_atomnames_matches_first_surface_pattern(tmp_path):
    surface_gro = tmp_path / "surface.gro"
    surface_itp = tmp_path / "surface.itp"
    surface_gro.write_text(
        "Surface\n"
        "    4\n"
        "    1SRF     P4    1   0.100   0.100   0.100\n"
        "    1SRF     C1    2   0.200   0.100   0.100\n"
        "    2SRF     P4    3   0.100   0.200   0.100\n"
        "    2SRF     C1    4   0.200   0.200   0.100\n"
        "   1.00000   1.00000   1.00000\n"
    )
    surface_itp.write_text(
        "[ moleculetype ]\n"
        "SRF 1\n\n"
        "[ atoms ]\n"
        "1 P4 1 SRF B1 1 0.0\n"
        "2 C1 1 SRF B2 2 0.0\n"
    )

    changed = pipeline._normalize_surface_itp_atomnames(surface_gro=surface_gro, surface_itp=surface_itp)

    assert changed is True
    text = surface_itp.read_text()
    assert "1 P4 1 SRF P4 1 0.0" in text
    assert "2 C1 1 SRF C1 2 0.0" in text


def test_sanitize_surface_itp_removes_legacy_marker_lines(tmp_path):
    surface_itp = tmp_path / "surface.itp"
    surface_itp.write_text(
        "[ moleculetype ]\n"
        "GRA 1\n\n"
        "[ atoms ]\n"
        "1 C1 1 GRA C 1 0.0\n"
        "~\n"
    )

    changed = pipeline._sanitize_surface_itp(surface_itp)

    assert changed is True
    assert "~" not in surface_itp.read_text()


def test_parser_rejects_ionize_without_solvate():
    parser = pipeline.build_parser()
    with pytest.raises(SystemExit):
        args = parser.parse_args([
            "--pdb", "1RJW",
            "--anchor", "1", "1",
            "--lx", "10",
            "--ly", "10",
            "--ionize",
        ])
        pipeline._validate_args(parser, args)


def test_parser_rejects_negative_salt_concentration():
    parser = pipeline.build_parser()
    with pytest.raises(SystemExit):
        args = parser.parse_args([
            "--pdb", "1RJW",
            "--anchor", "1", "1",
            "--lx", "10",
            "--ly", "10",
            "--solvate",
            "--salt-conc", "-0.1",
        ])
        pipeline._validate_args(parser, args)


def test_parser_rejects_non_positive_solvate_radius():
    parser = pipeline.build_parser()
    with pytest.raises(SystemExit):
        args = parser.parse_args([
            "--pdb", "1RJW",
            "--anchor", "1", "1",
            "--lx", "10",
            "--ly", "10",
            "--solvate",
            "--solvate-radius", "0.0",
        ])
        pipeline._validate_args(parser, args)


def test_parser_rejects_negative_solvate_surface_clearance():
    parser = pipeline.build_parser()
    with pytest.raises(SystemExit):
        args = parser.parse_args([
            "--pdb", "1RJW",
            "--anchor", "1", "1",
            "--lx", "10",
            "--ly", "10",
            "--solvate",
            "--solvate-surface-clearance", "-0.01",
        ])
        pipeline._validate_args(parser, args)


def test_parser_accepts_polarizable_water_for_dna():
    parser = pipeline.build_parser()
    args = parser.parse_args([
        "--dna",
        "--pdb", "4C64",
        "--anchor", "1", "1",
        "--lx", "10",
        "--ly", "10",
        "--solvate",
        "--polarizable-water",
    ])
    pipeline._validate_args(parser, args)
    assert args.polarizable_water is True


def test_parser_rejects_polarizable_water_without_dna():
    parser = pipeline.build_parser()
    with pytest.raises(SystemExit):
        args = parser.parse_args([
            "--pdb", "1RJW",
            "--anchor", "1", "1",
            "--lx", "10",
            "--ly", "10",
            "--solvate",
            "--polarizable-water",
        ])
        pipeline._validate_args(parser, args)


def test_write_ions_mdp_supports_polarizable_water(tmp_path):
    mdp = tmp_path / "ions.mdp"
    pipeline._write_ions_mdp(mdp, polarizable_water=True)
    text = mdp.read_text()

    assert "coulombtype              = Cut-off" in text
    assert "coulomb-modifier         = Potential-shift" in text
    assert "vdwtype                  = Cut-off" in text
    assert "vdw-modifier             = Force-switch" in text
    assert "constraints              = none" in text
    assert "epsilon-r                = 2.5" in text
    assert "verlet-buffer-tolerance = -1" in text


def test_write_ions_mdp_omits_epsilon_rf_for_standard_dna(tmp_path):
    mdp = tmp_path / "ions_dna.mdp"
    pipeline._write_ions_mdp(mdp, is_dna=True)
    text = mdp.read_text()

    assert "epsilon_r = 15" in text
    assert "epsilon_rf" not in text


def test_write_ions_mdp_keeps_epsilon_rf_for_non_dna(tmp_path):
    mdp = tmp_path / "ions_protein.mdp"
    pipeline._write_ions_mdp(mdp, is_dna=False)
    text = mdp.read_text()

    assert "epsilon_rf = 0" in text


def test_parser_accepts_dna_freeze_water_options():
    parser = pipeline.build_parser()
    args = parser.parse_args([
        "--dna",
        "--pdb", "4C64",
        "--anchor", "1", "1",
        "--lx", "10",
        "--ly", "10",
        "--solvate",
        "--freeze-water-fraction", "0.1",
        "--freeze-water-seed", "7",
    ])
    pipeline._validate_args(parser, args)
    assert args.freeze_water_fraction == 0.1
    assert args.freeze_water_seed == 7


def test_parser_rejects_freeze_water_without_dna():
    parser = pipeline.build_parser()
    with pytest.raises(SystemExit):
        args = parser.parse_args([
            "--pdb", "1RJW",
            "--anchor", "1", "1",
            "--lx", "10",
            "--ly", "10",
            "--solvate",
            "--freeze-water-fraction", "0.1",
        ])
        pipeline._validate_args(parser, args)


def test_run_genion_with_fallback_uses_second_attempt(monkeypatch, tmp_path):
    calls = {"n": 0}

    def fake_run_capture(cmd, cwd=None, stdin_text=None):
        calls["n"] += 1
        if calls["n"] == 1:
            return subprocess.CompletedProcess(cmd, 1, stdout="fail", stderr="bad group")
        return subprocess.CompletedProcess(cmd, 0, stdout="ok", stderr="")

    monkeypatch.setattr(pipeline, "_run_capture", fake_run_capture)

    pipeline._run_genion_with_fallback(
        gmx_bin="gmx",
        tpr_path=tmp_path / "ions.tpr",
        out_gro=tmp_path / "final.gro",
        top_path=tmp_path / "system_final.top",
        salt_conc=0.15,
    )
    assert calls["n"] == 2


def test_optional_solvation_requires_gromacs(monkeypatch, tmp_path):
    monkeypatch.setattr(pipeline, "_find_gmx_binary", lambda: None)
    args = Namespace(solvate=True, ionize=False, salt_conc=0.15, solvate_radius=0.21, solvate_surface_clearance=0.4)
    with pytest.raises(RuntimeError, match="no GROMACS binary"):
        pipeline._run_optional_solvation_ionization(args, tmp_path)


def test_sync_final_restrained_topology_copies_molecules_block(tmp_path):
    top_dir = tmp_path / "0_topology"
    top_dir.mkdir(parents=True)
    final_top = top_dir / "system_final.top"
    res_top = top_dir / "system_res.top"

    final_top.write_text(
        '#include "system_itp/martini_v3.0.0.itp"\n\n'
        "[ system ]\nX\n\n"
        "[ molecules ]\n"
        "Protein 1\n"
        "SRF 3\n"
        "W 100\n"
    )
    res_top.write_text(
        '#include "system_itp/martini_v3.0.0.itp"\n\n'
        "[ system ]\nX restrained\n\n"
        "[ molecules ]\n"
        "Protein 1\n"
        "SRF 1\n"
    )

    out = pipeline._sync_final_restrained_topology(top_dir, final_top)
    assert out is not None
    text = out.read_text()
    assert "SRF 3" in text
    assert "W 100" in text


def test_refresh_dna_thermostat_groups_rewrites_mdps_from_final_top(tmp_path):
    simdir = tmp_path / "Simulation"
    top_dir = simdir / "0_topology"
    mdp_dir = simdir / "1_mdp"
    itp_dir = top_dir / "system_itp"
    top_dir.mkdir(parents=True)
    mdp_dir.mkdir()
    itp_dir.mkdir()

    (top_dir / "system_final.top").write_text(
        "[ molecules ]\n"
        "Nucleic_A+Nucleic_B 1\n"
        "SRF 1\n"
        "W 20\n"
        "NA 2\n"
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
    (mdp_dir / "production_dna.mdp").write_text(
        "integrator = md\n"
        "rvdw                     = 1.1\n"
        "tcoupl                   = v-rescale\n"
        "tc-grps                  = System\n"
        "tau-t                    = 1.0\n"
        "ref-t                    = 300\n"
    )

    pipeline._refresh_dna_thermostat_groups(simdir=simdir, top_path=top_dir / "system_final.top")

    text = (mdp_dir / "production_dna.mdp").read_text()
    assert "define                   = -DPOSRES -DPOSRES_DNA" in text
    assert "tc-grps                  = SRF DNA W IONS" in text
    assert "tau-t                    = 0.3 0.3 0.3 0.3" in text
    assert "ref-t                    = 300 300 300 300" in text


def test_remove_waters_near_surface_updates_top(tmp_path):
    gro = tmp_path / "solvated_system.gro"
    top = tmp_path / "system_final.top"
    gro.write_text(
        "Test\n"
        "    6\n"
        "    1PRO    B1    1   0.100   0.100   1.000\n"
        "    2SRF    C1    2   0.500   0.500   0.300\n"
        "    3W       W    3   0.200   0.200   0.350\n"
        "    3W       W    4   0.200   0.200   0.350\n"
        "    4W       W    5   0.800   0.800   1.000\n"
        "    4W       W    6   0.800   0.800   1.000\n"
        "   2.00000   2.00000   2.00000\n"
    )
    top.write_text(
        "[ molecules ]\n"
        "Protein 1\n"
        "SRF 1\n"
        "W 2\n"
    )

    pipeline._remove_waters_below_surface(
        gro_path=gro,
        top_path=top,
        surface_resname="SRF",
        clearance_nm=0.21,
        water_resnames={"W"},
        water_molname="W",
    )
    t = top.read_text()
    assert "W               1" in t


def test_convert_standard_waters_to_polarizable_updates_gro_and_top(tmp_path):
    gro = tmp_path / "solvated_system.gro"
    top = tmp_path / "system_final.top"
    gro.write_text(
        "Test\n"
        "    3\n"
        "    1PRO    B1    1   0.100   0.100   1.000\n"
        "    2W       W    2   0.200   0.200   0.350\n"
        "    3W       W    3   0.800   0.800   1.000\n"
        "   2.00000   2.00000   2.00000\n"
    )
    top.write_text(
        "[ molecules ]\n"
        "Protein 1\n"
        "W 2\n"
    )

    converted = pipeline._convert_standard_waters_to_polarizable(gro, top, {"W"})

    assert converted == 2
    text = gro.read_text()
    assert "    7\n" in text
    assert "PW       W" in text
    assert "PW      WP" in text
    assert "PW      WM" in text
    top_text = top.read_text()
    assert "PW 2" in top_text
    assert "\nW 2\n" not in top_text


def test_rebuild_merged_index_writes_clean_groups_and_keeps_anchor_last(tmp_path):
    top_dir = tmp_path / "0_topology"
    (top_dir / "system_itp").mkdir(parents=True)
    gro = tmp_path / "2_system" / "system_final.gro"
    gro.parent.mkdir(parents=True)
    gro.write_text(
        "Dummy\n"
        "    4\n"
        "    1PRO     BB    1   0.100   0.100   0.100\n"
        "    2SRF     C1    2   0.200   0.200   0.200\n"
        "    3W       W     3   0.300   0.300   0.300\n"
        "    4NA     NA     4   0.400   0.400   0.400\n"
        "   1.00000   1.00000   1.00000\n"
    )
    (top_dir / "system_itp" / "surface.itp").write_text("[ moleculetype ]\nSRF 1\n\n[ atoms ]\n1 C1 1 SRF C1 1 0.0\n")

    index = top_dir / "index.ndx"
    index.write_text("[ LINKER ]\n2\n[ Anchor_1 ]\n1 2 3\n")

    out = pipeline._rebuild_merged_index(gmx_bin="gmx", gro_path=gro, top_dir=top_dir)
    assert out == index

    merged = index.read_text()
    assert "[ system ]" in merged
    assert "[ SRF ]" in merged
    assert "[ W ]" in merged
    assert "[ IONS ]" in merged
    assert "[ Protein ]" in merged
    assert "[ LINKER ]" in merged
    assert "[ Anchor_1 ]" in merged
    assert merged.index("[ LINKER ]") < merged.index("[ Anchor_1 ]")


def test_rebuild_merged_index_adds_extra_resname_groups_like_mol1(tmp_path):
    top_dir = tmp_path / "0_topology"
    (top_dir / "system_itp").mkdir(parents=True)
    gro = tmp_path / "2_system" / "system_final.gro"
    gro.parent.mkdir(parents=True)
    fmt = "{resid:5d}{resname:<5}{atomname:>5}{atomid:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n"
    gro.write_text(
        "Dummy\n"
        "    5\n"
        + fmt.format(resid=1, resname="DA", atomname="BB1", atomid=1, x=0.100, y=0.100, z=0.100)
        + fmt.format(resid=2, resname="MOL1", atomname="C01", atomid=2, x=0.200, y=0.200, z=0.200)
        + fmt.format(resid=2, resname="MOL1", atomname="C02", atomid=3, x=0.300, y=0.300, z=0.300)
        + fmt.format(resid=3, resname="SRF", atomname="C1", atomid=4, x=0.400, y=0.400, z=0.400)
        + fmt.format(resid=4, resname="W", atomname="W", atomid=5, x=0.500, y=0.500, z=0.500)
        + "   1.00000   1.00000   1.00000\n"
    )
    (top_dir / "system_itp" / "surface.itp").write_text("[ moleculetype ]\nSRF 1\n\n[ atoms ]\n1 C1 1 SRF C1 1 0.0\n")

    pipeline._rebuild_merged_index(gmx_bin="gmx", gro_path=gro, top_dir=top_dir)
    merged = (top_dir / "index.ndx").read_text()

    assert "[ DNA ]" in merged
    assert "[ MOL1 ]" in merged
    assert "[ W ]" in merged
    assert "[ SRF ]" in merged


def test_rebuild_merged_index_replaces_case_conflicting_auto_group(tmp_path):
    top_dir = tmp_path / "0_topology"
    (top_dir / "system_itp").mkdir(parents=True)
    gro = tmp_path / "2_system" / "system_final.gro"
    gro.parent.mkdir(parents=True)
    gro.write_text(
        "Dummy\n"
        "    2\n"
        "    1PRO     BB    1   0.100   0.100   0.100\n"
        "    2SRF     C1    2   0.200   0.200   0.200\n"
        "   1.00000   1.00000   1.00000\n"
    )
    (top_dir / "system_itp" / "surface.itp").write_text("[ moleculetype ]\nSRF 1\n\n[ atoms ]\n1 C1 1 SRF C1 1 0.0\n")

    index = top_dir / "index.ndx"
    index.write_text("[ System ]\n1 2\n[ system ]\n1 2\n[ Anchor_1 ]\n1 2\n")

    pipeline._rebuild_merged_index(gmx_bin="gmx", gro_path=gro, top_dir=top_dir)
    merged = index.read_text()

    assert "[ system ]" in merged
    assert "[ System ]" not in merged
    assert "[ Anchor_1 ]" in merged


def test_rebuild_merged_index_deduplicates_existing_custom_group_names(tmp_path):
    top_dir = tmp_path / "0_topology"
    (top_dir / "system_itp").mkdir(parents=True)
    gro = tmp_path / "2_system" / "system_final.gro"
    gro.parent.mkdir(parents=True)
    fmt = "{resid:5d}{resname:<5}{atomname:>5}{atomid:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n"
    gro.write_text(
        "Dummy\n"
        "    3\n"
        + fmt.format(resid=1, resname="DA", atomname="BB1", atomid=1, x=0.100, y=0.100, z=0.100)
        + fmt.format(resid=2, resname="MOL1", atomname="C01", atomid=2, x=0.200, y=0.200, z=0.200)
        + fmt.format(resid=2, resname="MOL1", atomname="C02", atomid=3, x=0.300, y=0.300, z=0.300)
        + "   1.00000   1.00000   1.00000\n"
    )
    (top_dir / "system_itp" / "surface.itp").write_text("[ moleculetype ]\nSRF 1\n\n[ atoms ]\n1 C1 1 SRF C1 1 0.0\n")

    index = top_dir / "index.ndx"
    index.write_text("[ MOL1 ]\n2 3\n[ MOL1 ]\n2 3\n[ Anchor_1 ]\n1 2\n")

    pipeline._rebuild_merged_index(gmx_bin="gmx", gro_path=gro, top_dir=top_dir)
    merged = index.read_text()

    assert merged.count("[ MOL1 ]") == 1
    assert merged.count("[ Anchor_1 ]") == 1


def test_restore_non_solvent_from_reference_keeps_solute_labels(tmp_path):
    ref = tmp_path / "solvated.gro"
    tgt = tmp_path / "final.gro"
    ref.write_text(
        "Ref\n"
        "    5\n"
        "    1PRO     BB    1   1.000   1.000   1.000\n"
        "    1PRO    SC1    2   1.100   1.100   1.100\n"
        "    2SRF     C1    3   0.500   0.500   0.100\n"
        "    3W       W     4   2.000   2.000   2.000\n"
        "    4W       W     5   2.500   2.500   2.500\n"
        "   3.00000   3.00000   3.00000\n"
    )
    tgt.write_text(
        "Final\n"
        "    5\n"
        "    1PRO     BB    1   9.000   9.000   9.000\n"
        "    1PRO    CA     2   9.100   9.100   9.100\n"
        "    2SRF     C1    3   9.200   9.200   9.200\n"
        "    3W       W     4   2.000   2.000   2.000\n"
        "    5NA      NA    5   2.500   2.500   2.500\n"
        "   3.00000   3.00000   3.00000\n"
    )

    ok = pipeline._restore_non_solvent_from_reference(
        reference_gro=ref,
        target_gro=tgt,
        water_resnames={"W"},
        ion_resnames={"NA", "CL"},
    )
    assert ok is True
    out = tgt.read_text()
    assert "PRO     BB" in out
    assert "PRO    SC1" in out
    assert "PRO    CA" not in out
    assert "NA      NA" in out
