from pathlib import Path

from martinisurf.utils.freeze_water import apply_freeze_water_fraction


def test_apply_freeze_water_fraction_updates_gro_and_top(tmp_path):
    gro = tmp_path / "final_system.gro"
    top = tmp_path / "system_final.top"
    alias = tmp_path / "system_final_alias.gro"

    atoms = [
        (1, "PRO", "B1", 1, 0.100, 0.100, 1.000),
        (2, "W", "W", 2, 0.200, 0.200, 1.100),
        (3, "W", "W", 3, 0.300, 0.300, 1.200),
        (4, "W", "W", 4, 0.400, 0.400, 1.300),
        (5, "W", "W", 5, 0.500, 0.500, 1.400),
        (6, "W", "W", 6, 0.600, 0.600, 1.500),
    ]
    atom_lines = [
        f"{resid:5d}{resname:<5}{atomname:>5}{atomid:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n"
        for resid, resname, atomname, atomid, x, y, z in atoms
    ]
    gro.write_text(
        "Test\n"
        "    6\n"
        + "".join(atom_lines)
        + "   2.00000   2.00000   2.00000\n"
    )
    top.write_text(
        "[ molecules ]\n"
        "Protein 1\n"
        "W 5\n"
    )

    n_before, n_w, n_wf = apply_freeze_water_fraction(
        top_path=top,
        gro_path=gro,
        fraction=0.2,
        seed=42,
        source_resname="W",
        target_resname="WF",
        alias_gro_path=alias,
    )

    assert n_before == 5
    assert n_w == 4
    assert n_wf == 1
    assert alias.exists()
    body = gro.read_text().splitlines()[2:-1]
    water_resnames = [ln[5:10].strip() for ln in body if ln[5:10].strip() in {"W", "WF"}]
    assert water_resnames == ["W", "W", "W", "W", "WF"]
    top_text = top.read_text()
    assert "W               4" in top_text
    assert "WF              1" in top_text
    assert top_text.index("W               4") < top_text.index("WF              1")
    frozen_lines = [ln for ln in body if ln[5:10].strip() == "WF"]
    assert len(frozen_lines) == 1
    assert frozen_lines[0][10:15].strip() == "WF"


def test_apply_freeze_water_fraction_spreads_wf_and_keeps_them_after_w(tmp_path):
    gro = tmp_path / "final_system.gro"
    top = tmp_path / "system_final.top"

    atoms = [
        (1, "W", "W", 1, 0.100, 0.100, 0.100),
        (2, "W", "W", 2, 0.900, 0.100, 0.100),
        (3, "W", "W", 3, 0.100, 0.900, 0.100),
        (4, "W", "W", 4, 0.900, 0.900, 0.100),
        (5, "W", "W", 5, 0.500, 0.500, 0.900),
    ]
    atom_lines = [
        f"{resid:5d}{resname:<5}{atomname:>5}{atomid:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n"
        for resid, resname, atomname, atomid, x, y, z in atoms
    ]
    gro.write_text(
        "Test\n"
        "    5\n"
        + "".join(atom_lines)
        + "   1.00000   1.00000   1.00000\n"
    )
    top.write_text(
        "[ molecules ]\n"
        "W 5\n"
    )

    _, n_w, n_wf = apply_freeze_water_fraction(
        top_path=top,
        gro_path=gro,
        fraction=0.4,
        seed=12345,
        source_resname="W",
        target_resname="WF",
    )

    assert (n_w, n_wf) == (3, 2)
    body = gro.read_text().splitlines()[2:-1]
    water_resnames = [ln[5:10].strip() for ln in body]
    frozen_resids = [int(ln[0:5]) for ln in body if ln[5:10].strip() == "WF"]

    assert water_resnames == ["W", "W", "W", "WF", "WF"]
    assert frozen_resids == [2, 5]
    assert frozen_resids != [4, 5]
    frozen_atomnames = [ln[10:15].strip() for ln in body if ln[5:10].strip() == "WF"]
    assert frozen_atomnames == ["WF", "WF"]


def test_apply_freeze_water_fraction_noop_with_zero_fraction(tmp_path):
    gro = tmp_path / "final_system.gro"
    top = tmp_path / "system_final.top"
    gro.write_text("T\n    0\n   1.0   1.0   1.0\n")
    top.write_text("[ molecules ]\n")

    out = apply_freeze_water_fraction(top_path=top, gro_path=gro, fraction=0.0)
    assert out == (0, 0, 0)
