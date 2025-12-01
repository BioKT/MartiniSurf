#!/usr/bin/env python3
"""
Utilities for downloading and cleaning PDB structures from RCSB
or AlphaFold UniProt database for MartiniSurf.
"""

from pathlib import Path
import urllib.request


RCSB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
AF_URL   = "https://alphafold.ebi.ac.uk/files/{uniprot_id}.pdb"


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


# ======================================================================
# Fetch from AlphaFold (UniProt ID)
# ======================================================================
def fetch_alphafold_pdb(uniprot_id: str, outdir: Path) -> Path:
    uniprot_id = uniprot_id.upper()
    outpath = outdir / f"{uniprot_id}_AF.pdb"
    url = AF_URL.format(uniprot_id=uniprot_id)

    print(f"⬇️  Downloading AlphaFold model for UniProt {uniprot_id}...")
    try:
        urllib.request.urlretrieve(url, outpath)
    except Exception as e:
        raise ValueError(f"Failed to download AlphaFold model for {uniprot_id}.\n{url}\n{e}")
    print(f"✔ Downloaded AlphaFold → {outpath}")
    return outpath


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
