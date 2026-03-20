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
    assert args.solvate is True
    assert args.ionize is False
    assert args.salt_conc == 0.20
    assert args.solvate_radius == 0.21
    assert args.solvate_surface_clearance == 0.4


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

    assert "coulombtype              = reaction-field" in text
    assert "vdw_type                 = cutoff" in text
    assert "vdw-modifier             = Force-switch" in text
    assert "constraints              = none" in text
    assert "epsilon_r                = 2.5" in text
    assert "epsilon_rf               = 0" in text


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


def test_rebuild_merged_index_appends_custom_groups(monkeypatch, tmp_path):
    top_dir = tmp_path / "0_topology"
    top_dir.mkdir(parents=True)
    gro = tmp_path / "2_system" / "system_final.gro"
    gro.parent.mkdir(parents=True)
    gro.write_text("dummy\n")

    index = top_dir / "index.ndx"
    index.write_text("[ Anchor_1 ]\n1 2 3\n")

    calls = {"n": 0}

    def fake_run_with_check(cmd, cwd=None, stdin_text=None):
        calls["n"] += 1
        assert cmd[0] == "gmx"
        assert cmd[1] == "make_ndx"
        assert stdin_text == "q\n"
        out_path = Path(cmd[cmd.index("-o") + 1])
        out_path.write_text("[ System ]\n1 2 3 4\n")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(pipeline, "_run_with_check", fake_run_with_check)

    out = pipeline._rebuild_merged_index(gmx_bin="gmx", gro_path=gro, top_dir=top_dir)
    assert out == index
    assert calls["n"] == 1

    merged = index.read_text()
    assert "[ System ]" in merged
    assert "[ Anchor_1 ]" in merged
    assert merged.index("[ System ]") < merged.index("[ Anchor_1 ]")


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
