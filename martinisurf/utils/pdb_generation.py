#!/usr/bin/env python3
"""
Utilities for downloading and cleaning PDB structures from RCSB
or AlphaFold UniProt database for MartiniSurf.
"""

import json
from pathlib import Path
import urllib.request
from urllib.error import URLError


RCSB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
AF_API_URL = "https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
AF_PDB_URL_PATTERNS = (
    "https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb",
    "https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v3.pdb",
    "https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v2.pdb",
    "https://alphafold.ebi.ac.uk/files/{uniprot_id}.pdb",
)


# ======================================================================
# Fetch from RCSB
# ======================================================================
def fetch_pdb(pdb_id: str, outdir: Path) -> Path:
    pdb_id = pdb_id.upper()
    outpath = outdir / f"{pdb_id}.pdb"
    url = RCSB_URL.format(pdb_id=pdb_id)

    print(f"⬇️  Downloading PDB {pdb_id} from RCSB...")
    try:
        urllib.request.urlretrieve(url, outpath)
    except Exception as e:
        raise ValueError(f"Failed to download RCSB ID {pdb_id}.\n{url}\n{e}")
    print(f"✔ Downloaded → {outpath}")
    return outpath


def _download_url(url: str, outpath: Path) -> Path:
    urllib.request.urlretrieve(url, outpath)
    return outpath


def _resolve_alphafold_pdb_url(uniprot_id: str) -> str | None:
    api_url = AF_API_URL.format(uniprot_id=uniprot_id)
    try:
        with urllib.request.urlopen(api_url) as response:
            payload = json.loads(response.read().decode("utf-8"))
    except Exception:
        return None

    if isinstance(payload, list):
        records = payload
    elif isinstance(payload, dict):
        records = [payload]
    else:
        return None

    for record in records:
        if not isinstance(record, dict):
            continue
        pdb_url = record.get("pdbUrl")
        if isinstance(pdb_url, str) and pdb_url:
            return pdb_url
        cif_url = record.get("cifUrl")
        if isinstance(cif_url, str) and cif_url:
            # When only mmCIF is reported, infer the sibling PDB file path.
            stem = cif_url.rsplit("/", 1)[-1]
            if stem.endswith(".cif"):
                return cif_url[:-4] + ".pdb"
    return None


# ======================================================================
# Fetch from AlphaFold (UniProt ID)
# ======================================================================
def fetch_alphafold_pdb(uniprot_id: str, outdir: Path) -> Path:
    uniprot_id = uniprot_id.upper()
    outpath = outdir / f"{uniprot_id}_AF.pdb"

    print(f"⬇️  Downloading AlphaFold model for UniProt {uniprot_id}...")
    tried_urls = []

    api_resolved_url = _resolve_alphafold_pdb_url(uniprot_id)
    if api_resolved_url:
        tried_urls.append(api_resolved_url)
        try:
            _download_url(api_resolved_url, outpath)
            print(f"✔ Downloaded AlphaFold → {outpath}")
            return outpath
        except Exception:
            pass

    for pattern in AF_PDB_URL_PATTERNS:
        url = pattern.format(uniprot_id=uniprot_id)
        if url in tried_urls:
            continue
        tried_urls.append(url)
        try:
            _download_url(url, outpath)
            print(f"✔ Downloaded AlphaFold → {outpath}")
            return outpath
        except Exception:
            continue

    raise ValueError(
        f"Failed to download AlphaFold model for {uniprot_id}.\n"
        f"Tried URLs:\n - " + "\n - ".join(tried_urls)
    )


# ======================================================================
# Cleaning
# ======================================================================
def simple_clean_pdb(infile: Path, outfile: Path, chain: str | None = None) -> Path:
    with open(infile) as fin, open(outfile, "w") as fout:
        for line in fin:
            if not line.startswith("ATOM"):
                continue
            if chain is not None:
                if len(line) < 22 or line[21] != chain:
                    continue
            fout.write(line)
    print(f"🧹 Cleaned PDB → {outfile}")
    return outfile


# ======================================================================
# Resolve local / RCSB / AlphaFold
# ======================================================================
def resolve_pdb_input(pdb_input: str, workdir: Path) -> Path:
    candidate = Path(pdb_input)
    if candidate.exists():
        return candidate.resolve()

    if len(pdb_input) == 4 and pdb_input.isalnum():
        return fetch_pdb(pdb_input, outdir=workdir)

    if len(pdb_input) == 6 and pdb_input.isalnum():
        return fetch_alphafold_pdb(pdb_input, outdir=workdir)

    raise ValueError(
        f"Invalid --pdb '{pdb_input}'. Must be:\n"
        " • local file\n"
        " • 4-letter RCSB PDB ID\n"
        " • 6-letter UniProt ID (AlphaFold)\n"
    )


# ======================================================================
# Main loader
# ======================================================================
def load_clean_pdb(pdb_input: str, workdir: Path, chain: str | None = None) -> Path:
    raw_pdb = resolve_pdb_input(pdb_input, workdir)

    cleaned = workdir / "2_system/cleaned_input.pdb"

    return simple_clean_pdb(raw_pdb, cleaned, chain=chain).resolve()
