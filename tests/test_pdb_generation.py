import json
from pathlib import Path

import pytest

from martinisurf.utils import pdb_generation


class _FakeResponse:
    def __init__(self, payload: str):
        self._payload = payload.encode("utf-8")

    def read(self):
        return self._payload

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False


def test_resolve_alphafold_pdb_url_prefers_api_pdb_url(monkeypatch):
    payload = json.dumps(
        [{"pdbUrl": "https://alphafold.ebi.ac.uk/files/AF-P69905-F1-model_v4.pdb"}]
    )

    def fake_urlopen(url):
        assert url == pdb_generation.AF_API_URL.format(uniprot_id="P69905")
        return _FakeResponse(payload)

    monkeypatch.setattr(pdb_generation.urllib.request, "urlopen", fake_urlopen)

    assert (
        pdb_generation._resolve_alphafold_pdb_url("P69905")
        == "https://alphafold.ebi.ac.uk/files/AF-P69905-F1-model_v4.pdb"
    )


def test_fetch_alphafold_pdb_falls_back_to_known_file_patterns(monkeypatch, tmp_path):
    monkeypatch.setattr(pdb_generation, "_resolve_alphafold_pdb_url", lambda _id: None)
    seen = []

    def fake_download(url, outpath):
        seen.append(url)
        if url.endswith("AF-P69905-F1-model_v4.pdb"):
            Path(outpath).write_text("ATOM\n")
            return outpath
        raise RuntimeError("missing")

    monkeypatch.setattr(pdb_generation, "_download_url", fake_download)

    outpath = pdb_generation.fetch_alphafold_pdb("P69905", tmp_path)

    assert outpath == tmp_path / "P69905_AF.pdb"
    assert outpath.read_text() == "ATOM\n"
    assert seen[0].endswith("AF-P69905-F1-model_v4.pdb")


def test_simple_clean_pdb_balances_merge_group_and_preserves_residue_ids(tmp_path):
    infile = tmp_path / "input.pdb"
    outfile = tmp_path / "output.pdb"
    infile.write_text(
        "ATOM      1  CA  ALA A   3      0.000   0.000   0.000  1.00 20.00           C\n"
        "ATOM      2  CA  GLY A   4      0.000   0.000   0.000  1.00 20.00           C\n"
        "ATOM      3  CA  SER A   5      0.000   0.000   0.000  1.00 20.00           C\n"
        "ATOM      4  CA  GLY B   4      0.000   0.000   0.000  1.00 20.00           C\n"
        "ATOM      5  CA  SER B   5      0.000   0.000   0.000  1.00 20.00           C\n"
        "ATOM      6  CA  THR B   6      0.000   0.000   0.000  1.00 20.00           C\n"
    )

    pdb_generation.simple_clean_pdb(
        infile,
        outfile,
        merge_groups=["A,B"],
        balance_merged_chains=True,
    )

    residues = [line[21:27] for line in outfile.read_text().splitlines()]
    assert residues == ["A   4 ", "A   5 ", "B   4 ", "B   5 "]


def test_simple_clean_pdb_fails_when_merge_group_remains_misaligned(tmp_path):
    infile = tmp_path / "input.pdb"
    outfile = tmp_path / "output.pdb"
    infile.write_text(
        "ATOM      1  CA  ALA A   3      0.000   0.000   0.000  1.00 20.00           C\n"
        "ATOM      2  CA  GLY A   4      0.000   0.000   0.000  1.00 20.00           C\n"
        "ATOM      3  CA  SER A   5      0.000   0.000   0.000  1.00 20.00           C\n"
        "ATOM      4  CA  GLY B   4      0.000   0.000   0.000  1.00 20.00           C\n"
        "ATOM      5  CA  SER B   5      0.000   0.000   0.000  1.00 20.00           C\n"
        "ATOM      6  CA  THR B   6      0.000   0.000   0.000  1.00 20.00           C\n"
    )

    with pytest.raises(ValueError, match="must expose the same ordered residue identities"):
        pdb_generation.simple_clean_pdb(
            infile,
            outfile,
            merge_groups=["A,B"],
            balance_merged_chains=False,
        )


def test_load_clean_pdb_converts_local_cif_before_cleaning(monkeypatch, tmp_path):
    cif = tmp_path / "input.cif"
    cif.write_text("data_test\n")

    def fake_convert_cif_to_pdb(cif_path, pdb_path):
        assert cif_path == cif
        pdb_path.write_text(
            "ATOM      1  CA  ALA A   1      0.000   0.000   0.000  1.00 20.00           C\n"
            "ATOM      2  CA  GLY A   2      0.100   0.000   0.000  1.00 20.00           C\n"
        )
        return pdb_path

    monkeypatch.setattr(pdb_generation, "convert_cif_to_pdb", fake_convert_cif_to_pdb)

    cleaned = pdb_generation.load_clean_pdb(str(cif), tmp_path)

    assert cleaned == (tmp_path / "2_system" / "cleaned_input.pdb").resolve()
    assert "ALA A   1" in cleaned.read_text()
