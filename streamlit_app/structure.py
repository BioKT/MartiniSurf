from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass
class ChainSummary:
    chain: str
    residues: int
    atoms: int


@dataclass
class StructureSummary:
    source: str
    format: str
    molecule_type: str
    chains: list[ChainSummary]
    residue_count: int
    atom_count: int
    ligands: list[str]
    notes: list[str]


AA3 = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS",
    "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "HSD", "HSE", "HSP", "MSE",
}
NUCLEIC = {"DA", "DT", "DG", "DC", "A", "U", "G", "C", "T"}
COMMON_SOLVENT = {"HOH", "WAT", "SOL", "NA", "CL", "K", "MG", "CA", "ZN"}


def summarize_structure(path: Path | None, identifier: str = "") -> StructureSummary:
    if path and path.exists():
        suffix = path.suffix.lower()
        if suffix == ".pdb":
            return _summarize_pdb(path)
        if suffix in {".cif", ".mmcif"}:
            return _summarize_cif(path)
        return StructureSummary(path.name, suffix.lstrip("."), "Unknown", [], 0, 0, [], ["Unsupported preview format."])

    clean_id = identifier.strip()
    if clean_id:
        return StructureSummary(
            clean_id,
            "remote",
            "Deferred",
            [],
            0,
            0,
            [],
            ["Remote identifiers are inspected during the MartiniSurf build step."],
        )

    return StructureSummary("No structure selected", "none", "Not configured", [], 0, 0, [], [])


def _summarize_pdb(path: Path) -> StructureSummary:
    residues: set[tuple[str, str, str]] = set()
    chain_atoms: dict[str, int] = {}
    chain_residues: dict[str, set[tuple[str, str]]] = {}
    residue_names: set[str] = set()
    ligands: set[str] = set()
    atom_count = 0

    for line in path.read_text(errors="replace").splitlines():
        record = line[:6].strip()
        if record not in {"ATOM", "HETATM"}:
            continue
        atom_count += 1
        resname = line[17:20].strip() or "UNK"
        chain = line[21].strip() or "_"
        resid = line[22:26].strip()
        icode = line[26].strip()
        residues.add((chain, resid, icode))
        residue_names.add(resname)
        chain_atoms[chain] = chain_atoms.get(chain, 0) + 1
        chain_residues.setdefault(chain, set()).add((resid, icode))
        if record == "HETATM" and resname not in COMMON_SOLVENT:
            ligands.add(resname)

    molecule_type = _guess_molecule_type(residue_names)
    chains = [
        ChainSummary(chain=chain, residues=len(chain_residues[chain]), atoms=chain_atoms[chain])
        for chain in sorted(chain_atoms)
    ]
    notes = [] if atom_count else ["No ATOM/HETATM records found."]
    return StructureSummary(path.name, "pdb", molecule_type, chains, len(residues), atom_count, sorted(ligands), notes)


def _summarize_cif(path: Path) -> StructureSummary:
    chains: dict[str, set[str]] = {}
    atoms = 0
    residue_names: set[str] = set()
    ligands: set[str] = set()

    for raw in path.read_text(errors="replace").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("ATOM ") or line.startswith("HETATM "):
            parts = line.split()
            if len(parts) < 9:
                continue
            atoms += 1
            record = parts[0]
            resname = parts[5] if len(parts) > 5 else "UNK"
            chain = parts[6] if len(parts) > 6 else "_"
            resid = parts[8] if len(parts) > 8 else str(atoms)
            residue_names.add(resname)
            chains.setdefault(chain, set()).add(resid)
            if record == "HETATM" and resname not in COMMON_SOLVENT:
                ligands.add(resname)

    chain_summary = [ChainSummary(chain=chain, residues=len(resids), atoms=0) for chain, resids in sorted(chains.items())]
    notes = [] if atoms else ["mmCIF preview is limited; full validation happens in MartiniSurf."]
    return StructureSummary(
        path.name,
        "mmcif",
        _guess_molecule_type(residue_names),
        chain_summary,
        sum(len(resids) for resids in chains.values()),
        atoms,
        sorted(ligands),
        notes,
    )


def _guess_molecule_type(residue_names: set[str]) -> str:
    if not residue_names:
        return "Unknown"
    protein_hits = len(residue_names & AA3)
    nucleic_hits = len(residue_names & NUCLEIC)
    if protein_hits and nucleic_hits:
        return "Protein with nucleic-acid-like residues"
    if nucleic_hits > protein_hits:
        return "Non-protein-like"
    if protein_hits:
        return "Protein-like"
    return "Unknown"
