from pathlib import Path

import pytest

from martinisurf.utils.use_python2 import find_python2


def _write_fake_python(path: Path, major: int) -> None:
    path.write_text(
        "#!/usr/bin/env bash\n"
        "if [ \"$1\" = \"-c\" ]; then\n"
        f"  echo {major}\n"
        "  exit 0\n"
        "fi\n"
        "exit 1\n"
    )
    path.chmod(0o755)


def test_find_python2_prefers_real_python2_candidate(tmp_path):
    py3 = tmp_path / "py3fake"
    py2 = tmp_path / "py2fake"
    _write_fake_python(py3, major=3)
    _write_fake_python(py2, major=2)

    selected = find_python2(candidates=[str(py3), str(py2)])
    assert selected == str(py2)


def test_find_python2_raises_when_no_python2(tmp_path):
    py3 = tmp_path / "py3fake"
    _write_fake_python(py3, major=3)

    with pytest.raises(RuntimeError) as exc:
        find_python2(candidates=[str(py3)])

    assert "requires python2.7" in str(exc.value)
