from __future__ import annotations

import pytest

from streamlit_app.short_md import (
    ShortMDConfig,
    ShortMDStage,
    extract_mdrun_performance,
    patch_mdp_text,
    preview_rows,
    safe_nsteps,
    validate_stage_order,
)


def test_safe_nsteps_matches_colab_formula() -> None:
    assert safe_nsteps(time_ns=0.02, dt_ps=0.001) == 20000
    assert safe_nsteps(time_ns=0.1, dt_ps=0.01) == 10000


def test_safe_nsteps_rejects_non_positive_values() -> None:
    with pytest.raises(ValueError, match="Time step"):
        safe_nsteps(time_ns=0.1, dt_ps=0.0)
    with pytest.raises(ValueError, match="Time"):
        safe_nsteps(time_ns=0.0, dt_ps=0.01)


def test_patch_mdp_text_updates_existing_and_adds_missing_keys() -> None:
    text = """
integrator               = md
dt                       = 0.020000
nsteps                   = 10
""".lstrip()

    patched = patch_mdp_text(text, dt_ps=0.01, nsteps=5000, nstxoutc=1000)

    assert "dt                       = 0.010000" in patched
    assert "nsteps                   = 5000" in patched
    assert "nstxout-compressed       = 1000" in patched


def test_validate_stage_order_matches_colab_dependencies() -> None:
    stages = [ShortMDStage("production", True, 0.01, 0.1)]
    assert validate_stage_order(stages) == ["Deposition and Production require NPT to be enabled."]

    stages = [
        ShortMDStage("npt", True, 0.005, 0.02),
        ShortMDStage("production", True, 0.01, 0.1),
    ]
    assert validate_stage_order(stages) == ["NPT requires NVT to be enabled."]


def test_preview_rows_uses_xtc_stride_from_colab() -> None:
    config = ShortMDConfig(
        xtc_write_every_ps=10.0,
        stages=(ShortMDStage("production", True, 0.01, 0.1),),
    )

    assert preview_rows(config) == [
        {
            "Stage": "PRODUCTION",
            "dt (ps)": 0.01,
            "time (ns)": 0.1,
            "nsteps": 10000,
            "saved frame dt (ps)": 10.0,
        }
    ]


def test_extract_mdrun_performance() -> None:
    text = "Performance:      123.456 ns/day,      0.194 hours/ns"
    assert extract_mdrun_performance(text) == (123.456, 0.194)
    assert extract_mdrun_performance("no performance line") == (None, None)
