#!/usr/bin/env python3
"""
Utilities for downloading and cleaning PDB structures from the RCSB database
for use with MartiniSurf and Martinize2.

Implements the minimal cleaning recommended by the official Martinize2
tutorial:
    • Only keep lines beginning with 'ATOM'
    • Optionally keep a single protein chain (e.g. chain A)
"""

from pathlib import Path
import urllib.request


# Base URL for RCSB PDB download
RCSB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"


# ======================================================================
# Download PDB from RCSB
# ======================================================================
def fetch_pdb(pdb_id: str, outdir: Path) -> Path:
    """
    Download a PDB structure from the RCSB server.

    Parameters
    ----------
    pdb_id : str
        A 4-character PDB ID (case-insensitive)
    outdir : Path
        Directory where the file should be written

    Returns
    -------
    Path
        Full path to downloaded PDB file
    """
    pdb_id = pdb_id.lower()
    outpath = outdir / f"{pdb_id}.pdb"
    url = RCSB_URL.format(pdb_id=pdb_id)

    print(f"⬇️  Downloading PDB {pdb_id.upper()} from RCSB...")

    try:
        urllib.request.urlretrieve(url, outpath)
    except Exception as e:
        raise ValueError(
            f"Failed to download PDB ID {pdb_id.upper()} from RCSB.\nURL: {url}\n{e}"
        )

    print(f"✔ Downloaded → {outpath}")
    return outpath


# ======================================================================
# Minimal cleaning (Martinize2 tutorial)
# ======================================================================
def simple_clean_pdb(infile: Path, outfile: Path, chain: str | None = None) -> Path:
    """
    Minimal PDB cleaning following Martinize2 recommendations:
      • Keep ONLY ATOM records.
      • Optionally filter by chain (chain A, B, etc.).

    Parameters
    ----------
    infile : Path
        Raw PDB file
    outfile : Path
        Output cleaned file
    chain : str or None
        Chain identifier to keep ('A', 'B', etc.)
        If None, keep all chains.

    Returns
    -------
    Path
        Path to cleaned PDB file
    """
    with open(infile, "r") as fin, open(outfile, "w") as fout:
        for line in fin:
            if not line.startswith("ATOM"):
                continue

            if chain is not None:
                # PDB chain column is index 21
                if len(line) < 22 or line[21] != chain:
                    continue

            fout.write(line)

    print(f"🧹 Cleaned PDB → {outfile}")
    return outfile


# ======================================================================
# Resolve PDB input (local file OR PDB ID)
# ======================================================================
def resolve_pdb_input(pdb_input: str, workdir: Path) -> Path:
    """
    Resolve user input for --pdb:
        • If it's a local file → return absolute path.
        • If it's a 4-letter PDB ID → download from RCSB.
        • Otherwise → error.

    Parameters
    ----------
    pdb_input : str
        User-supplied text from --pdb
    workdir : Path
        Directory where downloaded files should be saved

    Returns
    -------
    Path
        Path to local PDB file
    """
    # local file
    candidate = Path(pdb_input)
    if candidate.exists():
        return candidate.resolve()

    # RCSB PDB ID
    if len(pdb_input) == 4 and pdb_input.isalnum():
        return fetch_pdb(pdb_input, outdir=workdir)

    # otherwise error
    raise ValueError(
        f"Invalid PDB input '{pdb_input}'. Must be:\n"
        "  • Local PDB file path, OR\n"
        "  • Valid 4-letter PDB ID (e.g. 1CRN, 4JJR)"
    )


# ======================================================================
# Final combined loader used by the pipeline
# ======================================================================
def load_clean_pdb(pdb_input: str, workdir: Path, chain: str | None = None) -> Path:
    """
    Resolve and clean a PDB for MartiniSurf:
        • Accepts local file or PDB ID
        • Downloads if needed
        • Cleans using simple_clean_pdb

    Parameters
    ----------a
    pdb_input : str
        The string passed via --pdb
    workdir : Path
        Where cleaned PDB should be written
    chain : str or None
        Select specific chain (optional)

    Returns
    -------
    Path
        Absolute path to cleaned PDB ready for Martinize2
    """
    raw_pdb = resolve_pdb_input(pdb_input, workdir)
    cleaned = workdir / "2_system/cleaned_input.pdb"
    return simple_clean_pdb(raw_pdb, cleaned, chain=chain).resolve()
