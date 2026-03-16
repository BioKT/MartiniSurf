import json
from pathlib import Path

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
