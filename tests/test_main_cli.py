import sys

import martinisurf
import martinisurf.__main__ as cli


def test_explicit_build_subcommand_forwards_flags_without_build_token(monkeypatch):
    captured = {}

    class DummyModule:
        @staticmethod
        def main():
            captured["argv"] = sys.argv[:]

    monkeypatch.setattr(cli, "_print_compact_package_help", lambda: None)
    monkeypatch.setattr(cli, "_print_full_package_help", lambda: None)
    monkeypatch.setitem(sys.modules, "martinisurf.pipeline", DummyModule)
    monkeypatch.setattr(martinisurf, "pipeline", DummyModule, raising=False)
    monkeypatch.setattr(sys, "argv", ["martinisurf", "build", "--pdb", "1RJW"])

    cli.main()

    assert captured["argv"] == ["martinisurf build", "--pdb", "1RJW"]
