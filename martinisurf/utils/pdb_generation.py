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
def _is_cif_path(path: Path) -> bool:
    return path.suffix.lower() in {".cif", ".mmcif"}


def convert_cif_to_pdb(cif_path: Path, pdb_path: Path) -> Path:
    """Convert an mmCIF/PDBx structure to PDB for the existing cleaning path."""
    try:
        import mdtraj as md
    except ImportError as exc:
        raise RuntimeError(
            "Reading .cif/.mmcif inputs requires MDTraj, which is listed as a MartiniSurf dependency."
        ) from exc

    print(f"🔄 Converting mmCIF → PDB (MDTraj): {cif_path} → {pdb_path}")
    pdb_path.parent.mkdir(parents=True, exist_ok=True)
    traj = md.load(str(cif_path))
    traj.save_pdb(str(pdb_path))
    print(f"✔ Converted mmCIF → {pdb_path}")
    return pdb_path


def _parse_merge_groups(raw_groups: list[str] | None, chain_ids: list[str]) -> list[list[str]]:
    if not raw_groups:
        return []

    parsed_groups: list[list[str]] = []
    for raw in raw_groups:
        item = str(raw).strip()
        if not item:
            continue
        if item.lower() == "all":
            if len(chain_ids) > 1:
                parsed_groups.append(list(chain_ids))
            continue

        members = [part.strip() for part in item.split(",") if part.strip()]
        members = [member for member in members if member in chain_ids]
        if len(members) > 1:
            parsed_groups.append(members)
    return parsed_groups


def _residue_key_from_pdb_line(line: str) -> tuple[int, str, str]:
    resid = int(line[22:26])
    icode = line[26].strip() if len(line) > 26 else ""
    resname = line[17:20].strip() if len(line) >= 20 else ""
    return resid, icode, resname


def _build_chain_residue_order(
    atom_lines: list[str],
) -> dict[str, list[tuple[int, str, str]]]:
    chain_residue_order: dict[str, list[tuple[int, str, str]]] = {}
    seen_per_chain: dict[str, set[tuple[int, str, str]]] = {}

    for line in atom_lines:
        if len(line) < 26:
            continue
        chain_id = (line[21].strip() or "_")
        residue_key = _residue_key_from_pdb_line(line)
        chain_residue_order.setdefault(chain_id, [])
        seen_per_chain.setdefault(chain_id, set())
        if residue_key not in seen_per_chain[chain_id]:
            chain_residue_order[chain_id].append(residue_key)
            seen_per_chain[chain_id].add(residue_key)

    return chain_residue_order


def validate_merged_chain_residue_alignment(
    atom_lines: list[str],
    merge_groups: list[str] | None,
) -> dict[str, list[tuple[int, str, str]]]:
    chain_residue_order = _build_chain_residue_order(atom_lines)
    parsed_groups = _parse_merge_groups(merge_groups, list(chain_residue_order))
    validated_groups: dict[str, list[tuple[int, str, str]]] = {}
    for group in parsed_groups:
        reference_chain = group[0]
        reference_keys = list(chain_residue_order.get(reference_chain, []))
        for chain_id in group[1:]:
            chain_keys = list(chain_residue_order.get(chain_id, []))
            if chain_keys != reference_keys:
                reference_preview = ", ".join(
                    f"{resname}{resid}{icode}" for resid, icode, resname in reference_keys[:8]
                ) or "none"
                chain_preview = ", ".join(
                    f"{resname}{resid}{icode}" for resid, icode, resname in chain_keys[:8]
                ) or "none"
                raise ValueError(
                    "Merged chains must expose the same ordered residue identities before any feature merge. "
                    f"Group {','.join(group)} is misaligned between chains {reference_chain} and {chain_id}. "
                    f"Reference residues: [{reference_preview}] | {chain_id}: [{chain_preview}]. "
                    "Re-enable automatic balancing or fix the input numbering so merged chains match exactly."
                )
        validated_groups[",".join(group)] = reference_keys
    return validated_groups


def _build_balanced_residue_selection(
    atom_lines: list[str],
    merge_groups: list[str] | None,
) -> dict[str, set[tuple[int, str, str]]]:
    chain_residue_order = _build_chain_residue_order(atom_lines)
    parsed_groups = _parse_merge_groups(merge_groups, list(chain_residue_order))
    if not parsed_groups:
        return {}

    allowed_by_chain: dict[str, set[tuple[int, str, str]]] = {}
    for group in parsed_groups:
        ordered_counts = {chain_id: len(chain_residue_order.get(chain_id, [])) for chain_id in group}
        ref_chain = min(group, key=lambda chain_id: ordered_counts.get(chain_id, 0))
        residue_sets = {chain_id: set(chain_residue_order.get(chain_id, [])) for chain_id in group}
        common_keys = [
            key for key in chain_residue_order.get(ref_chain, [])
            if all(key in residue_sets[other_chain] for other_chain in group)
        ]
        allowed = set(common_keys)
        for chain_id in group:
            allowed_by_chain[chain_id] = allowed
        print(
            "🧹 Balanced merged chains "
            f"{','.join(group)} using reference {ref_chain}: "
            f"{ordered_counts.get(ref_chain, 0)} -> {len(common_keys)} shared residues"
        )

    return allowed_by_chain


def simple_clean_pdb(
    infile: Path,
    outfile: Path,
    chain: str | None = None,
    merge_groups: list[str] | None = None,
    balance_merged_chains: bool = False,
    validate_merged_alignment: bool = True,
) -> Path:
    atom_lines: list[str] = []
    with open(infile) as fin:
        for line in fin:
            if not line.startswith("ATOM"):
                continue
            if chain is not None and (len(line) < 22 or line[21] != chain):
                continue
            atom_lines.append(line)

    allowed_by_chain = (
        _build_balanced_residue_selection(atom_lines, merge_groups)
        if balance_merged_chains else {}
    )
    filtered_lines: list[str] = []

    for line in atom_lines:
        chain_id = (line[21].strip() or "_")
        if chain_id in allowed_by_chain:
            if _residue_key_from_pdb_line(line) not in allowed_by_chain[chain_id]:
                continue
        filtered_lines.append(line)

    if validate_merged_alignment:
        validate_merged_chain_residue_alignment(filtered_lines, merge_groups)

    with open(outfile, "w") as fout:
        for line in filtered_lines:
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

    system_dir = workdir / "2_system"

    if len(pdb_input) == 4 and pdb_input.isalnum():
        return fetch_pdb(pdb_input, outdir=system_dir)

    if len(pdb_input) == 6 and pdb_input.isalnum():
        return fetch_alphafold_pdb(pdb_input, outdir=system_dir)

    raise ValueError(
        f"Invalid --pdb '{pdb_input}'. Must be:\n"
        " • local file\n"
        " • 4-letter RCSB PDB ID\n"
        " • 6-letter UniProt ID (AlphaFold)\n"
    )


# ======================================================================
# Main loader
# ======================================================================
def load_clean_pdb(
    pdb_input: str,
    workdir: Path,
    chain: str | None = None,
    merge_groups: list[str] | None = None,
    balance_merged_chains: bool = False,
    validate_merged_alignment: bool = True,
) -> Path:
    raw_structure = resolve_pdb_input(pdb_input, workdir)
    raw_pdb = raw_structure
    if _is_cif_path(raw_structure):
        cif_pdb = workdir / "2_system" / f"{raw_structure.stem}_from_cif.pdb"
        cif_pdb.parent.mkdir(parents=True, exist_ok=True)
        raw_pdb = convert_cif_to_pdb(
            raw_structure,
            cif_pdb,
        )

    cleaned = workdir / "2_system/cleaned_input.pdb"

    return simple_clean_pdb(
        raw_pdb,
        cleaned,
        chain=chain,
        merge_groups=merge_groups,
        balance_merged_chains=balance_merged_chains,
        validate_merged_alignment=validate_merged_alignment,
    ).resolve()
