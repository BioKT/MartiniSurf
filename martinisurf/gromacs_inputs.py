#!/usr/bin/env python3
"""
GoMartini System Builder with Anchoring Restraints

Original stable version
+ Safe linker support
"""

from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path
from typing import Dict, List, Sequence

import MDAnalysis as mda
import martinisurf


DNA_STANDARD_FF = "martini_v2.1-dna.itp"
DNA_POLARIZABLE_FF = "martini_v2.1P-dna.itp"
STANDARD_WATER_TEMPLATE = "water.gro"
POLARIZABLE_WATER_TEMPLATE = "polarize-water.gro"


# ======================================================================
# UTILITIES
# ======================================================================

def ensure_dir(path: str | Path) -> None:
    Path(path).mkdir(parents=True, exist_ok=True)


def write_list(values: List[int], fh, chunk: int = 15) -> None:
    for i in range(0, len(values), chunk):
        fh.write(" ".join(str(v) for v in values[i: i + chunk]) + "\n")


def _select_dna_forcefield_name(itp_dir: Path, polarizable_water: bool = False) -> str | None:
    preferred = DNA_POLARIZABLE_FF if polarizable_water else DNA_STANDARD_FF
    fallback = DNA_STANDARD_FF if polarizable_water else DNA_POLARIZABLE_FF

    if (itp_dir / preferred).exists():
        return preferred
    if polarizable_water:
        raise FileNotFoundError(
            f"Polarizable-water mode requires {preferred} in {itp_dir}."
        )
    if (itp_dir / fallback).exists():
        return fallback
    return None


def _polarizable_water_electrostatics_block() -> list[str]:
    return [
        "; OPTIONS FOR ELECTROSTATICS AND VDW =",
        "; Martini 2 DNA polarizable-water setup modernized for GROMACS 2024 =",
        "; Legacy Coulomb Shift is approximated with Cut-off + Potential-shift =",
        "coulombtype              = Cut-off",
        "coulomb-modifier         = Potential-shift",
        "rcoulomb                 = 1.2",
        "; Dielectric constant =",
        "epsilon-r                = 2.5",
        "; Legacy VdW Shift is approximated with Cut-off + Force-switch =",
        "vdwtype                  = Cut-off",
        "vdw-modifier             = Force-switch",
        "; cut-off lengths        =",
        "rvdw-switch              = 0.9",
        "rvdw                     = 1.2",
        "; Apply long range dispersion corrections for Energy and Pressure =",
        "DispCorr                 = No",
    ]


def _polarizable_water_bond_block() -> list[str]:
    return [
        "; OPTIONS FOR BONDS     =",
        "constraints              = none",
        "; Type of constraint algorithm =",
        "constraint_algorithm     = Lincs",
        "; Relative tolerance of shake =",
        "shake_tol                = 0.0001",
        "; Highest order in the expansion of the constraint coupling matrix =",
        "lincs_order              = 4",
        "; Lincs will write a warning to the stderr if in one step a bond =",
        "; rotates over more degrees than =",
        "lincs_warnangle          = 90",
    ]


def _rewrite_mdp_for_polarizable_water(mdp_path: Path) -> None:
    lines = mdp_path.read_text().splitlines()

    electrostatic_keys = {
        "coulombtype",
        "coulomb-modifier",
        "rcoulomb_switch",
        "rcoulomb",
        "epsilon_r",
        "epsilon-r",
        "epsilon_rf",
        "vdw_type",
        "vdwtype",
        "vdw-modifier",
        "rvdw_switch",
        "rvdw-switch",
        "rvdw",
        "dispcorr",
    }
    bond_keys = {
        "constraints",
        "constraint_algorithm",
        "shake_tol",
        "lincs_order",
        "lincs_warnangle",
    }

    kept: list[str] = []
    insert_after = -1
    for raw in lines:
        stripped = raw.strip()
        if stripped:
            key = stripped.split("=", 1)[0].strip().lower() if "=" in stripped else ""
            if key in electrostatic_keys or key in bond_keys:
                continue
            if (
                key == "verlet-buffer-tolerance"
                or key == "rlist"
                or stripped.lower().startswith("; nblist cut-off")
            ):
                insert_after = len(kept)
        kept.append(raw)

    if not any(line.strip().lower().startswith("verlet-buffer-tolerance") for line in kept):
        if insert_after >= 0:
            kept.insert(insert_after + 1, "verlet-buffer-tolerance  = -1")
        else:
            kept.insert(0, "verlet-buffer-tolerance  = -1")

    electro_block = _polarizable_water_electrostatics_block()
    if insert_after < 0:
        insert_at = len(kept)
    else:
        insert_at = insert_after + 1
        while insert_at < len(kept) and not kept[insert_at].strip():
            insert_at += 1

    out = kept[:insert_at]
    if out and out[-1].strip():
        out.append("")
    out.extend(electro_block)

    tail = kept[insert_at:]
    if tail and tail[0].strip():
        out.append("")
    out.extend(tail)

    while out and not out[-1].strip():
        out.pop()
    out.append("")
    out.extend(_polarizable_water_bond_block())
    mdp_path.write_text("\n".join(out).rstrip() + "\n")


def parse_anchor_groups(anchor_args: List[List[str]] | None) -> Dict[int, List[int]]:
    groups: Dict[int, List[int]] = {}
    if not anchor_args:
        return groups
    for group in anchor_args:
        gid = int(group[0])
        residues = list(map(int, group[1:]))
        groups[gid] = residues
    return dict(sorted(groups.items()))


def build_anchor_atoms(
    universe: mda.Universe,
    anchor_groups: Dict[int, List[int]],
    exclude_atom_ids: set[int],
) -> Dict[int, List[int]]:
    anchor_atoms: Dict[int, List[int]] = {}
    for gid, residues in anchor_groups.items():
        atom_list = [
            int(a.index + 1)
            for a in universe.atoms
            if a.resnum in residues and int(a.index + 1) not in exclude_atom_ids
        ]
        anchor_atoms[gid] = atom_list
    return anchor_atoms


def write_top_files(
    topo_dir: Path,
    dst_itp_dir: Path,
    moltype: str,
    mol_itp_name: str,
    anchor_itp_name: str,
    is_dna: bool,
    use_linker: bool,
    go_model: bool = False,
    polarizable_water: bool = False,
    linker_itp_name: str = "linker.itp",
    restrained_linker_itp_name: str | None = None,
    linker_moltype: str | None = None,
    linker_count: int = 0,
    cofactor_itp_name: str | None = None,
    cofactor_moltype: str | None = None,
    cofactor_count: int = 0,
    substrate_itp_name: str | None = None,
    substrate_moltype: str | None = None,
    substrate_count: int = 0,
    surface_moltype: str = "SRF",
    surface_count: int = 1,
    intermolecular_itp_name: str | None = None,
) -> None:
    """Create simple master topology files expected by legacy workflows/tests."""
    if is_dna:
        # DNA FF files are mutually exclusive (both define [ defaults ]).
        dna_ff = _select_dna_forcefield_name(
            dst_itp_dir,
            polarizable_water=polarizable_water,
        )
        forcefield_itps = [dna_ff, "martini_v2.0_ions.itp"] if dna_ff else ["martini_v2.0_ions.itp"]
    else:
        forcefield_itps = [
            "martini_v3.0.0.itp",
            "martini_v3.0.0_solvents_v1.itp",
            "martini_v3.0.0_ions_v1.itp",
        ]

    present_forcefields = [itp for itp in forcefield_itps if itp and (dst_itp_dir / itp).exists()]

    def _include_block(molecule_itp_name: str, linker_include_name: str | None = None) -> str:
        lines = []
        if is_dna:
            lines.append("#define RUBBER_BANDS")
            lines.append("")
        elif go_model:
            lines.append("#define GO_VIRT")
            lines.append("")
        lines.extend(f'#include "system_itp/{name}"' for name in present_forcefields)
        lines.append(f'#include "system_itp/{molecule_itp_name}"')
        lines.append('#include "system_itp/surface.itp"')
        if use_linker:
            chosen_linker_itp = linker_include_name or linker_itp_name
            lines.append(f'#include "system_itp/{chosen_linker_itp}"')
        if cofactor_itp_name and cofactor_count > 0:
            lines.append(f'#include "system_itp/{cofactor_itp_name}"')
        if substrate_itp_name and substrate_count > 0:
            lines.append(f'#include "system_itp/{substrate_itp_name}"')
        if intermolecular_itp_name:
            lines.append(f'#include "system_itp/{intermolecular_itp_name}"')
        return "\n".join(lines)

    system_top = topo_dir / "system.top"
    system_res_top = topo_dir / "system_res.top"
    system_anchor_top = topo_dir / "system_anchor.top"
    linker_line = ""
    if use_linker and linker_count > 0:
        linker_name = linker_moltype or Path(linker_itp_name).stem
        linker_line = f"{linker_name} {linker_count}\n"
    cofactor_line = ""
    if cofactor_count > 0:
        cofactor_name = cofactor_moltype or (Path(cofactor_itp_name).stem if cofactor_itp_name else "COF")
        cofactor_line = f"{cofactor_name} {cofactor_count}\n"
    substrate_line = ""
    if substrate_count > 0:
        substrate_name = substrate_moltype or (Path(substrate_itp_name).stem if substrate_itp_name else "SUB")
        substrate_line = f"{substrate_name} {substrate_count}\n"

    with open(system_top, "w") as fh:
        fh.write(
            _include_block(mol_itp_name, linker_itp_name)
            + "\n\n[ system ]\nMartiniSurf system\n\n[ molecules ]\n"
            + f"{moltype} 1\n{cofactor_line}{linker_line}{surface_moltype} {surface_count}\n{substrate_line}"
        )

    with open(system_res_top, "w") as fh:
        fh.write(
            _include_block(anchor_itp_name, restrained_linker_itp_name)
            + "\n\n[ system ]\nMartiniSurf restrained system\n\n[ molecules ]\n"
            + f"{moltype} 1\n{cofactor_line}{linker_line}{surface_moltype} {surface_count}\n{substrate_line}"
        )

    # Compatibility alias for legacy workflows/scripts that still expect system_anchor.top.
    shutil.copy(system_res_top, system_anchor_top)


def _read_itp_moleculetype(itp_path: Path) -> str | None:
    if not itp_path.exists():
        return None

    lines = itp_path.read_text().splitlines()
    in_moleculetype = False
    for raw in lines:
        line = raw.strip()
        if not line or line.startswith(";"):
            continue
        if line.lower().startswith("[") and "moleculetype" in line.lower():
            in_moleculetype = True
            continue
        if in_moleculetype:
            if line.startswith("["):
                break
            token = line.split()[0]
            return token
    return None


def _read_itp_moleculetype_set(itp_path: Path) -> set[str]:
    moltypes: set[str] = set()
    if not itp_path.exists():
        return moltypes

    lines = itp_path.read_text().splitlines()
    in_moleculetype = False
    for raw in lines:
        line = raw.strip()
        if not line or line.startswith(";"):
            continue
        if line.lower().startswith("[") and "moleculetype" in line.lower():
            in_moleculetype = True
            continue
        if in_moleculetype:
            if line.startswith("["):
                in_moleculetype = False
                continue
            token = line.split()[0]
            moltypes.add(token)
            in_moleculetype = False
    return moltypes


def _read_itp_first_atoms_resname(itp_path: Path) -> str | None:
    if not itp_path.exists():
        return None

    in_atoms = False
    for raw in itp_path.read_text().splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line:
            continue
        if line.startswith("["):
            section = line.strip("[]").strip().lower()
            in_atoms = section == "atoms"
            continue
        if not in_atoms:
            continue
        parts = line.split()
        if len(parts) >= 5:
            return parts[3]
    return None


def _write_itp_with_moleculetype(src_itp: Path, dst_itp: Path, new_moltype: str) -> None:
    lines = src_itp.read_text().splitlines()
    out: list[str] = []
    in_moleculetype = False
    replaced = False
    for raw in lines:
        line = raw.strip()
        if line.lower().startswith("[") and "moleculetype" in line.lower():
            in_moleculetype = True
            out.append(raw)
            continue
        if in_moleculetype:
            if (not line) or line.startswith(";"):
                out.append(raw)
                continue
            if line.startswith("["):
                in_moleculetype = False
                out.append(raw)
                continue
            # Replace first data line in [ moleculetype ]
            tokens = raw.split()
            if tokens:
                nrexcl = tokens[1] if len(tokens) > 1 else "1"
                out.append(f"{new_moltype} {nrexcl}")
                replaced = True
                in_moleculetype = False
                continue
        out.append(raw)

    if not replaced:
        dst_itp.write_text(src_itp.read_text())
    else:
        dst_itp.write_text("\n".join(out).rstrip() + "\n")


def _count_itp_atoms(itp_path: Path) -> int:
    if not itp_path.exists():
        return 0
    count = 0
    in_atoms = False
    for raw in itp_path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith("["):
            in_atoms = "atoms" in line.lower()
            continue
        if in_atoms:
            # Count only valid atom definition lines where first token is atom id.
            # This ignores stray markers (e.g., "~") that may appear in some legacy ITPs.
            parts = line.split()
            if not parts:
                continue
            if parts[0].lstrip("+-").isdigit():
                count += 1
    return count


def _rewrite_itp_with_posres(src_itp: Path, dst_itp: Path, posres_atom_ids: List[int]) -> None:
    """Copy an ITP replacing any existing [ position_restraints ] block with a generated one."""
    out_lines: list[str] = []
    in_posres = False

    for raw in src_itp.read_text().splitlines(keepends=True):
        stripped = raw.strip().lower()
        if stripped.startswith("[") and "position_restraints" in stripped:
            in_posres = True
            continue
        if in_posres:
            if raw.strip().startswith("["):
                in_posres = False
                out_lines.append(raw)
            continue
        out_lines.append(raw)

    if posres_atom_ids:
        if out_lines and not out_lines[-1].endswith("\n"):
            out_lines[-1] = out_lines[-1] + "\n"
        out_lines.append("\n[ position_restraints ]\n")
        out_lines.append("#ifdef POSRES\n")
        for atom_id in sorted(set(posres_atom_ids)):
            out_lines.append(f"{atom_id} 1 1000 1000 0\n")
        out_lines.append("#endif\n")

    dst_itp.write_text("".join(out_lines))


def _materialize_posres_fc_in_itp(itp_path: Path, default_fc: float = 1000.0) -> bool:
    """
    Replace POSRES_FC tokens in non-preprocessor lines with a numeric value.
    Uses #define POSRES_FC <value> when available; otherwise default_fc.
    """
    if not itp_path.exists():
        return False

    text = itp_path.read_text()
    if "POSRES_FC" not in text:
        return False

    define_match = re.search(
        r"^\s*#\s*define\s+POSRES_FC\s+([0-9eE+\-.]+)\s*$",
        text,
        flags=re.MULTILINE,
    )
    fc_value = define_match.group(1) if define_match else f"{float(default_fc):.1f}"

    out_lines: list[str] = []
    changed = False
    for raw in text.splitlines(keepends=True):
        stripped = raw.lstrip()
        if stripped.startswith("#"):
            out_lines.append(raw)
            continue
        new_raw = re.sub(r"\bPOSRES_FC\b", fc_value, raw)
        if new_raw != raw:
            changed = True
        out_lines.append(new_raw)

    if changed:
        itp_path.write_text("".join(out_lines))
    return changed


def _infer_linker_restrained_atom_ids(
    pull_anchor_atoms: Dict[int, List[int]],
    linker_atom_ids: set[int],
    linker_size: int | None,
) -> List[int]:
    """
    Infer linker-local atom indices to restrain from linker-tail pull groups.
    Falls back to the last linker bead when mapping is unavailable.
    """
    if not linker_size or linker_size <= 0:
        return []

    tail_global_ids = sorted(
        {
            int(atoms[0])
            for gid, atoms in pull_anchor_atoms.items()
            if gid % 3 == 0 and atoms
        }
    )
    if not tail_global_ids:
        return [linker_size]

    linker_ids_sorted = sorted(int(i) for i in linker_atom_ids)
    id_to_pos = {atom_id: idx for idx, atom_id in enumerate(linker_ids_sorted)}

    local_ids: set[int] = set()
    for tail_id in tail_global_ids:
        pos = id_to_pos.get(tail_id)
        if pos is None:
            continue
        local_ids.add((pos % linker_size) + 1)

    return sorted(local_ids) if local_ids else [linker_size]


def _pick_linker_neighbor(
    universe: mda.Universe,
    instance_atom_ids: List[int],
    linker_head_id: int,
) -> int | None:
    """Pick the linker bead immediately adjacent to the head (fallback: nearest bead)."""
    others = [idx for idx in instance_atom_ids if idx != linker_head_id]
    if not others:
        return None
    if len(others) == 1:
        return others[0]

    head_pos = universe.atoms[linker_head_id - 1].position
    return min(
        others,
        key=lambda idx: float(((universe.atoms[idx - 1].position - head_pos) ** 2).sum()),
    )


def _pick_preferred_dna_bead(
    universe: mda.Universe,
    selected_atom_ids: List[int],
    linker_head_id: int,
) -> int | None:
    """
    Select DNA bead by priority BB1 -> BB2 -> BB3.
    If multiple candidates exist, choose the closest to linker_head.
    """
    if not selected_atom_ids:
        return None

    head_pos = universe.atoms[linker_head_id - 1].position
    selected = [universe.atoms[idx - 1] for idx in selected_atom_ids]

    for bead_name in ("BB1", "BB2", "BB3"):
        candidates = [a for a in selected if str(a.name).strip().upper() == bead_name]
        if candidates:
            best = min(
                candidates,
                key=lambda a: float(((a.position - head_pos) ** 2).sum()),
            )
            return int(best.index + 1)

    return None


def _parse_itp_section_tokens(itp_path: Path, sections: set[str]) -> dict[str, list[list[str]]]:
    out: dict[str, list[list[str]]] = {name: [] for name in sections}
    current: str | None = None
    for raw in itp_path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(";"):
            continue
        if line.startswith("[") and line.endswith("]"):
            current_name = line[1:-1].strip().lower()
            current = current_name if current_name in sections else None
            continue
        if current is None:
            continue
        data = raw.split(";", 1)[0].strip()
        if not data:
            continue
        tokens = data.split()
        if tokens:
            out[current].append(tokens)
    return out


def _append_itp_entries(lines: list[str], section: str, entries: list[str]) -> list[str]:
    if not entries:
        return lines

    target = section.lower()
    start = -1
    for i, raw in enumerate(lines):
        stripped = raw.strip().lower()
        if stripped.startswith("[") and stripped.endswith("]"):
            if stripped[1:-1].strip() == target:
                start = i
                break

    if start < 0:
        out = list(lines)
        if out and out[-1].strip():
            out.append("")
        out.append(f"[ {section} ]")
        out.extend(entries)
        return out

    end = len(lines)
    for j in range(start + 1, len(lines)):
        stripped = lines[j].strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            end = j
            break

    return lines[:end] + entries + lines[end:]


def _merge_dna_linker_itp(
    dst_itp_path: Path,
    dna_itp_path: Path,
    linker_itp_path: Path,
    universe: mda.Universe,
    linker_pairs: List[dict[str, List[int] | int]],
    linker_resname: str,
    bond_length_nm: float,
    bond_force: float = 1250.0,
    angle_deg: float = 180.0,
    angle_force: float = 20.0,
) -> tuple[int, int]:
    """
    Merge linker atoms/interactions into DNA moleculetype and add DNA-linker bond/angle.
    Returns (n_bonds_added, n_angles_added) for DNA-linker coupling.
    """
    base_lines = dna_itp_path.read_text().splitlines()
    base_atom_count = _count_itp_atoms(dna_itp_path)
    if base_atom_count <= 0:
        raise ValueError(f"DNA topology has no [ atoms ] entries: {dna_itp_path}")

    parsed = _parse_itp_section_tokens(
        linker_itp_path,
        {"atoms", "bonds", "angles", "constraints", "dihedrals", "pairs"},
    )
    linker_atoms_tpl = parsed["atoms"]
    if not linker_atoms_tpl:
        raise ValueError(f"Linker topology has no [ atoms ] entries: {linker_itp_path}")

    dna_atoms_section = _parse_itp_section_tokens(dna_itp_path, {"atoms"})["atoms"]
    max_resnr = max((int(row[2]) for row in dna_atoms_section if len(row) > 2 and row[2].lstrip("+-").isdigit()), default=1)
    max_cgnr = max((int(row[5]) for row in dna_atoms_section if len(row) > 5 and row[5].lstrip("+-").isdigit()), default=base_atom_count)

    linker_global_ids = sorted(
        {
            int(gid)
            for pair in linker_pairs
            for gid in pair["instance_atoms"]  # type: ignore[index]
        }
    )
    global_to_local = {gid: base_atom_count + i + 1 for i, gid in enumerate(linker_global_ids)}

    atom_entries: list[str] = []
    bonds_entries: list[str] = []
    angles_entries: list[str] = []
    constraints_entries: list[str] = []
    dihedrals_entries: list[str] = []
    pairs_entries: list[str] = []
    dna_linker_bonds = 0
    dna_linker_angles = 0
    tail_posres_local_ids: set[int] = set()

    next_resnr = max_resnr + 1
    next_cgnr = max_cgnr + 1

    template_atoms = linker_atoms_tpl
    n_tpl = len(template_atoms)
    if n_tpl <= 0:
        raise ValueError("Invalid linker template atom count.")

    for pair in linker_pairs:
        instance_atoms = [int(v) for v in pair["instance_atoms"]]  # type: ignore[index]
        if len(instance_atoms) != n_tpl:
            raise ValueError(
                f"Linker atom count mismatch: topology has {n_tpl}, instance has {len(instance_atoms)}."
            )

        tpl_to_local = {
            int(template_atoms[i][0]): global_to_local[instance_atoms[i]]
            for i in range(n_tpl)
        }
        linker_tail_global = int(pair["linker_tail"])  # type: ignore[arg-type]
        if linker_tail_global in global_to_local:
            tail_posres_local_ids.add(global_to_local[linker_tail_global])

        for i, row in enumerate(template_atoms):
            local_id = tpl_to_local[int(row[0])]
            atom_type = row[1] if len(row) > 1 else "C1"
            atom_name = row[4] if len(row) > 4 else f"L{i+1}"
            charge = row[6] if len(row) > 6 else "0.0"
            mass = row[7] if len(row) > 7 else None
            entry = f"{local_id:6d} {atom_type:>6} {next_resnr:6d} {linker_resname:<6} {atom_name:>6} {next_cgnr:6d} {charge:>8}"
            if mass is not None:
                entry += f" {mass}"
            atom_entries.append(entry)
            next_cgnr += 1
        next_resnr += 1

        for row in parsed["bonds"]:
            if len(row) < 2:
                continue
            i_local = tpl_to_local[int(row[0])]
            j_local = tpl_to_local[int(row[1])]
            rest = " ".join(row[2:]) if len(row) > 2 else "1"
            bonds_entries.append(f"{i_local:6d} {j_local:6d} {rest}")

        for row in parsed["angles"]:
            if len(row) < 3:
                continue
            i_local = tpl_to_local[int(row[0])]
            j_local = tpl_to_local[int(row[1])]
            k_local = tpl_to_local[int(row[2])]
            rest = " ".join(row[3:]) if len(row) > 3 else "1 180.0 20.0"
            angles_entries.append(f"{i_local:6d} {j_local:6d} {k_local:6d} {rest}")

        for row in parsed["constraints"]:
            if len(row) < 2:
                continue
            i_local = tpl_to_local[int(row[0])]
            j_local = tpl_to_local[int(row[1])]
            rest = " ".join(row[2:]) if len(row) > 2 else "1 0.30"
            constraints_entries.append(f"{i_local:6d} {j_local:6d} {rest}")

        for row in parsed["dihedrals"]:
            if len(row) < 4:
                continue
            i_local = tpl_to_local[int(row[0])]
            j_local = tpl_to_local[int(row[1])]
            k_local = tpl_to_local[int(row[2])]
            l_local = tpl_to_local[int(row[3])]
            rest = " ".join(row[4:]) if len(row) > 4 else "1 0.0 0.0"
            dihedrals_entries.append(f"{i_local:6d} {j_local:6d} {k_local:6d} {l_local:6d} {rest}")

        for row in parsed["pairs"]:
            if len(row) < 2:
                continue
            i_local = tpl_to_local[int(row[0])]
            j_local = tpl_to_local[int(row[1])]
            rest = " ".join(row[2:]) if len(row) > 2 else "1"
            pairs_entries.append(f"{i_local:6d} {j_local:6d} {rest}")

        selected_atoms = [int(v) for v in pair["selected_atoms"]]  # type: ignore[index]
        linker_head = int(pair["linker_head"])  # type: ignore[arg-type]
        dna_bead = _pick_preferred_dna_bead(
            universe=universe,
            selected_atom_ids=selected_atoms,
            linker_head_id=linker_head,
        )
        if dna_bead is None:
            continue

        linker_head_local = global_to_local[linker_head]
        bonds_entries.append(
            f"{dna_bead:6d} {linker_head_local:6d} 1 {bond_length_nm:.3f} {bond_force:.1f}"
        )
        dna_linker_bonds += 1

        linker_neighbor = _pick_linker_neighbor(
            universe=universe,
            instance_atom_ids=instance_atoms,
            linker_head_id=linker_head,
        )
        if linker_neighbor is not None:
            linker_neighbor_local = global_to_local[linker_neighbor]
            angles_entries.append(
                f"{linker_neighbor_local:6d} {linker_head_local:6d} {dna_bead:6d} 1 {angle_deg:.1f} {angle_force:.1f}"
            )
            dna_linker_angles += 1

    merged = list(base_lines)
    merged = _append_itp_entries(merged, "atoms", atom_entries)
    merged = _append_itp_entries(merged, "bonds", bonds_entries)
    merged = _append_itp_entries(merged, "angles", angles_entries)
    merged = _append_itp_entries(merged, "constraints", constraints_entries)
    merged = _append_itp_entries(merged, "dihedrals", dihedrals_entries)
    merged = _append_itp_entries(merged, "pairs", pairs_entries)
    dst_itp_path.write_text("\n".join(merged).rstrip() + "\n")
    # Keep linker-tail restraint for the bead pulled against the surface.
    _rewrite_itp_with_posres(
        src_itp=dst_itp_path,
        dst_itp=dst_itp_path,
        posres_atom_ids=sorted(tail_posres_local_ids),
    )
    return dna_linker_bonds, dna_linker_angles


def _infer_surface_pull_group(text: str, is_dna: bool, surface_moltype: str | None = None) -> str:
    if surface_moltype:
        return str(surface_moltype).strip()
    m = re.search(r"^\s*pull-group2-name\s*=\s*([A-Za-z0-9_+-]+)", text, flags=re.MULTILINE)
    if m:
        return m.group(1)
    return "SRF"


def _strip_pull_section(text: str) -> str:
    keep = []
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped:
            keep.append(line)
            continue
        if stripped.startswith(";pull"):
            continue
        if stripped.startswith("pull"):
            continue
        keep.append(line)
    return "\n".join(keep).rstrip() + "\n"


def _build_pull_block(
    anchor_count: int,
    surface_group: str,
    k: float = 500.0,
    init_nm: float | None = None,
    include_init: bool = True,
) -> str:
    if anchor_count <= 0:
        return ""

    lines = [
        "",
        "pull                     = yes",
        f"pull-ngroups             = {anchor_count + 1}",
    ]

    for i in range(1, anchor_count + 1):
        lines.append(f"pull-group{i}-name         = Anchor_{i}")

    surface_gid = anchor_count + 1
    lines.append(f"pull-group{surface_gid}-name         = {surface_group}")
    lines.append("")
    lines.append(f"pull-ncoords             = {anchor_count}")

    for i in range(1, anchor_count + 1):
        lines.append(f"pull-coord{i}-geometry     = distance")
        lines.append(f"pull-coord{i}-groups       = {i} {surface_gid}")
        lines.append(f"pull-coord{i}-type         = umbrella")
        lines.append(f"pull-coord{i}-k            = {k:.1f}")
        lines.append(f"pull-coord{i}-rate         = 0")
        lines.append(f"pull-coord{i}-dim          = N N Y")
        _append_pull_start_or_init(
            lines=lines,
            coord_id=i,
            include_init=include_init,
            init_nm=init_nm,
        )
        lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def _append_pull_start_or_init(
    lines: list[str],
    coord_id: int,
    include_init: bool,
    init_nm: float | None,
) -> None:
    """Write either _start or -init, never both, for one pull coordinate."""
    has_init = include_init and init_nm is not None
    if not has_init:
        lines.append(f"pull-coord{coord_id}_start        = yes")
        return
    lines.append(f"pull-coord{coord_id}-init         = {init_nm:.3f}")


def _build_linker_pull_block(
    linker_count: int,
    surface_group: str,
    k: float = 500.0,
    init_prot_nm: float | None = None,
    init_surf_nm: float | None = None,
    include_init: bool = True,
    include_head_pull: bool = True,
) -> str:
    if linker_count <= 0:
        return ""

    surface_gid = 3 * linker_count + 1

    lines = [
        "",
        "pull                     = yes",
        f"pull-ngroups             = {surface_gid}",
    ]

    for gid in range(1, surface_gid):
        lines.append(f"pull-group{gid}-name         = Anchor_{gid}")
    lines.append(f"pull-group{surface_gid}-name         = {surface_group}")
    lines.append("")
    coord_total = 2 * linker_count if include_head_pull else linker_count
    lines.append(f"pull-ncoords             = {coord_total}")

    coord_id = 1
    for i in range(linker_count):
        selected_gid = 3 * i + 1
        linker_head_gid = 3 * i + 2
        linker_tail_gid = 3 * i + 3

        if include_head_pull:
            # coord: biomolecule anchor <-> linker head (near biomolecule)
            lines.append(f"pull-coord{coord_id}-geometry     = distance")
            lines.append(f"pull-coord{coord_id}-groups       = {selected_gid} {linker_head_gid}")
            lines.append(f"pull-coord{coord_id}-type         = umbrella")
            lines.append(f"pull-coord{coord_id}-k            = {k:.1f}")
            lines.append(f"pull-coord{coord_id}-rate         = 0")
            lines.append(f"pull-coord{coord_id}-dim          = Y Y Y")
            _append_pull_start_or_init(
                lines=lines,
                coord_id=coord_id,
                include_init=include_init,
                init_nm=init_prot_nm,
            )
            lines.append("")
            coord_id += 1

        # coord: linker tail (near surface) <-> surface
        lines.append(f"pull-coord{coord_id}-geometry     = distance")
        lines.append(f"pull-coord{coord_id}-groups       = {linker_tail_gid} {surface_gid}")
        lines.append(f"pull-coord{coord_id}-type         = umbrella")
        lines.append(f"pull-coord{coord_id}-k            = {k:.1f}")
        lines.append(f"pull-coord{coord_id}-rate         = 0")
        lines.append(f"pull-coord{coord_id}-dim          = N N Y")
        _append_pull_start_or_init(
            lines=lines,
            coord_id=coord_id,
            include_init=include_init,
            init_nm=init_surf_nm,
        )
        lines.append("")
        coord_id += 1

    return "\n".join(lines).rstrip() + "\n"


def write_custom_mdp(
    src: Path,
    dst: Path,
    anchor_count: int,
    is_dna: bool,
    surface_moltype: str | None = None,
    linker_pull: bool = False,
    linker_count: int = 1,
    linker_pull_init_prot_nm: float | None = None,
    linker_pull_init_surf_nm: float | None = None,
    anchor_pull_init_nm: float | None = None,
    rewrite_pull: bool = True,
    include_pull_init: bool = True,
    dna_linker_bonded: bool = False,
) -> None:
    text = src.read_text()
    if not rewrite_pull:
        dst.write_text(text)
        return
    has_pull = bool(re.search(r"^\s*pull\s*=", text, flags=re.MULTILINE))
    surface_group = _infer_surface_pull_group(text, is_dna=is_dna, surface_moltype=surface_moltype)
    clean = _strip_pull_section(text)

    if has_pull:
        if linker_pull:
            clean = clean.rstrip() + "\n" + _build_linker_pull_block(
                linker_count=linker_count,
                surface_group=surface_group,
                init_prot_nm=linker_pull_init_prot_nm,
                init_surf_nm=linker_pull_init_surf_nm,
                include_init=include_pull_init,
                include_head_pull=not dna_linker_bonded,
            )
        elif anchor_count > 0:
            clean = clean.rstrip() + "\n" + _build_pull_block(
                anchor_count,
                surface_group,
                init_nm=anchor_pull_init_nm,
                include_init=include_pull_init,
            )

    dst.write_text(clean)


def _validate_cli_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    linker_mode = _linker_mode_enabled(args)

    if linker_mode and args.linker_size is not None and args.linker_size <= 0:
        parser.error("--linker-size must be > 0.")
    if args.cofactor_count < 0:
        parser.error("--cofactor-count must be >= 0.")
    if args.cofactor_count > 0 and not args.cofactor_itp_name:
        parser.error("--cofactor-count requires --cofactor-itp-name.")
    if args.substrate_count < 0:
        parser.error("--substrate-count must be >= 0.")
    if args.substrate_count > 0 and not (args.substrate_itp_name or args.substrate_moltype):
        parser.error("--substrate-count requires --substrate-itp-name or --substrate-moltype.")
    if linker_mode and args.linker_pull_init_prot is not None and args.linker_pull_init_prot <= 0:
        parser.error("--linker-pull-init-prot must be > 0 when provided.")
    if linker_mode and args.linker_pull_init_surf is not None and args.linker_pull_init_surf <= 0:
        parser.error("--linker-pull-init-surf must be > 0 when provided.")
    if args.ads_mode and linker_mode:
        parser.error("--ads-mode is incompatible with linker mode.")


def _linker_mode_enabled(args: argparse.Namespace) -> bool:
    """Linker handling is only active when explicitly requested or using the legacy linker anchor mode."""
    return bool(args.use_linker or args.linker_resid is not None)


def _normalize_linker_args(args: argparse.Namespace) -> bool:
    """
    Ignore stray linker-related CLI values unless linker mode is actually active.
    Returns the effective linker mode.
    """
    linker_mode = _linker_mode_enabled(args)
    if linker_mode:
        return True

    args.linker_resname = None
    args.linker_size = None
    args.linker_pull_init_prot = None
    args.linker_pull_init_surf = None
    return False


# ======================================================================
# MAIN WORKFLOW
# ======================================================================

def main(argv: Sequence[str] | None = None) -> None:

    parser = argparse.ArgumentParser(
        description="Build GoMartini system with anchoring restraints.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("--moltype", required=False, help="Molecule name for protein systems.")
    parser.add_argument("--outdir", default="Simulation", help="Output root (contains 0_topology/1_mdp/2_system).")

    # 🔹 Classical anchors
    parser.add_argument(
        "--anchor",
        nargs="+",
        action="append",
        metavar=("GROUP", "RESID"),
        help="Anchor group definition: GROUP RESID [RESID ...]. Repeat for multiple groups.",
    )

    # 🔹 NEW — Linker mode
    parser.add_argument("--linker-resid", type=int, help="Legacy single anchor residue for linker mode.")
    parser.add_argument("--use-linker", action="store_true", help="Enable linker-aware pull/index generation.")
    parser.add_argument("--linker-resname", help="Residue name used by linker beads in immobilized_system.gro.")
    parser.add_argument("--linker-size", type=int, help="Number of beads per linker instance.")
    parser.add_argument("--linker-itp-name", default="linker.itp", help="Linker topology filename in system_itp.")
    parser.add_argument(
        "--linker-pull-init-prot",
        type=float,
        default=None,
        help="Initial pull distance (nm) for linker-head to biomolecule anchor coordinates.",
    )
    parser.add_argument(
        "--linker-pull-init-surf",
        type=float,
        default=None,
        help="Initial pull distance (nm) for linker-tail to surface coordinates.",
    )
    parser.add_argument("--go-model", action="store_true", help="Add GO_VIRT define to generated protein topologies.")
    parser.add_argument(
        "--ads-mode",
        action="store_true",
        help=(
            "Adsorption mode: skip anchor pull/posres generation and emit only "
            "minimization/NVT/NPT/production MDPs (no deposition)."
        ),
    )
    parser.add_argument("--cofactor-itp-name", help="Cofactor topology filename in system_itp.")
    parser.add_argument("--cofactor-count", type=int, default=0, help="Number of cofactor molecules in the system.")
    parser.add_argument("--substrate-itp-name", help="Substrate topology filename in system_itp.")
    parser.add_argument("--substrate-moltype", help="Substrate moleculetype/resname when defined by the Martini force field.")
    parser.add_argument("--substrate-count", type=int, default=0, help="Number of substrate molecules in the system.")
    parser.add_argument(
        "--polarizable-water",
        action="store_true",
        help="DNA-only: use polarizable water templates, PW solvent, and martini_v2.1P-dna.itp.",
    )

    args = parser.parse_args(argv) if argv else parser.parse_args()
    _validate_cli_args(parser, args)
    linker_mode = _normalize_linker_args(args)
    ads_mode = bool(args.ads_mode)

    # ===============================================================
    # Locate immobilized_system.gro
    # ===============================================================
    cwd = Path.cwd()
    candidates = [
        cwd / "2_system" / "immobilized_system.gro",
        cwd / "immobilized_system.gro",
        Path(args.outdir).resolve() / "2_system" / "immobilized_system.gro",
    ]

    input_gro = next((c for c in candidates if c.exists()), None)
    if input_gro is None:
        raise FileNotFoundError("❌ immobilized_system.gro not found.")

    output_root = input_gro.parent.parent.resolve()

    print(f"\n• Using input system: {input_gro}")
    print(f"• Output directory:   {output_root}\n")

    u = mda.Universe(str(input_gro))

    # ===============================================================
    # AUTO-DETECT DNA vs PROTEIN
    # ===============================================================
    dna_resnames = {"DA", "DT", "DC", "DG", "DAN", "THY", "GUA", "CYT"}
    is_dna = any(res in dna_resnames for res in u.residues.resnames)

    print("🔍 Molecule type detected:", "DNA" if is_dna else "PROTEIN/MARTINI3")
    if args.polarizable_water and not is_dna:
        raise ValueError("❌ --polarizable-water is currently supported only for DNA systems.")

    # ===============================================================
    # BUILD ANCHOR GROUPS
    # ===============================================================

    anchor_atoms: Dict[int, List[int]] = {}
    pull_anchor_atoms: Dict[int, List[int]] = {}
    linker_pull_enabled = False
    linker_atom_ids: set[int] = set()
    linker_pairs: List[dict[str, List[int] | int]] = []

    if linker_mode and args.linker_resname:
        linker_atom_ids = {
            int(a.index + 1)
            for a in u.atoms
            if str(a.resname).strip() == args.linker_resname
        }

    # 🔵 LINKER MODE
    if args.linker_resid is not None:

        print(f"\n🔗 Linker mode → using residue {args.linker_resid} as anchor")

        atom_list = [
            int(a.index + 1)
            for a in u.atoms
            if a.resnum == args.linker_resid and int(a.index + 1) not in linker_atom_ids
        ]

        if not atom_list:
            raise ValueError(f"❌ Residue {args.linker_resid} not found.")

        anchor_atoms[1] = atom_list
        pull_anchor_atoms[1] = atom_list
        print(f"  → Anchor_1: {len(atom_list)} atoms")

    # 🔵 CLASSICAL MODE
    else:
        anchor_groups = parse_anchor_groups(args.anchor)
        if not anchor_groups:
            print("ℹ No anchors provided. Anchor restraints will be skipped.")
        else:
            print("\n=== Anchor Groups Provided ===")
            for gid, res in anchor_groups.items():
                print(f"  Anchor_{gid}: {res}")
            print()

        anchor_atoms = build_anchor_atoms(u, anchor_groups, linker_atom_ids)
        pull_anchor_atoms.update(anchor_atoms)
        for gid, atom_list in anchor_atoms.items():
            print(f"  → Anchor_{gid}: {len(atom_list)} atoms")

    if linker_atom_ids and pull_anchor_atoms:
        selected_items = [(gid, atoms) for gid, atoms in pull_anchor_atoms.items() if atoms]
        linker_ids_sorted = sorted(linker_atom_ids)

        if args.linker_size and args.linker_size > 0:
            n_instances = len(linker_ids_sorted) // args.linker_size
            n_pairs = min(len(selected_items), n_instances)

            if n_pairs > 0:
                new_pull_groups: Dict[int, List[int]] = {}

                for i in range(n_pairs):
                    selected_atoms = selected_items[i][1]
                    instance = linker_ids_sorted[i * args.linker_size: (i + 1) * args.linker_size]
                    if not instance:
                        continue

                    selected_centroid = u.atoms[[idx - 1 for idx in selected_atoms]].positions.mean(axis=0)
                    linker_positions = {idx: u.atoms[idx - 1].position for idx in instance}

                    linker_head = min(
                        instance,
                        key=lambda idx: float(((linker_positions[idx] - selected_centroid) ** 2).sum()),
                    )
                    linker_tail = min(instance, key=lambda idx: float(linker_positions[idx][2]))

                    if linker_tail == linker_head and len(instance) > 1:
                        linker_tail = sorted(
                            (idx for idx in instance if idx != linker_head),
                            key=lambda idx: float(linker_positions[idx][2]),
                        )[0]

                    base = 3 * i + 1
                    new_pull_groups[base] = selected_atoms
                    new_pull_groups[base + 1] = [linker_head]
                    new_pull_groups[base + 2] = [linker_tail]
                    linker_pairs.append(
                        {
                            "selected_atoms": selected_atoms,
                            "instance_atoms": instance,
                            "linker_head": linker_head,
                            "linker_tail": linker_tail,
                        }
                    )

                    print(f"  → Anchor_{base + 1} (linker head): atom {linker_head}")
                    print(f"  → Anchor_{base + 2} (linker tail): atom {linker_tail}")

                if new_pull_groups:
                    pull_anchor_atoms = new_pull_groups
                    linker_pull_enabled = True

        else:
            # Backward-compatible fallback: treat all linker atoms as one linker instance.
            selected_atoms = selected_items[0][1]
            selected_centroid = u.atoms[[idx - 1 for idx in selected_atoms]].positions.mean(axis=0)
            linker_positions = {idx: u.atoms[idx - 1].position for idx in linker_ids_sorted}

            linker_head = min(
                linker_ids_sorted,
                key=lambda idx: float(((linker_positions[idx] - selected_centroid) ** 2).sum()),
            )
            linker_tail = min(linker_ids_sorted, key=lambda idx: float(linker_positions[idx][2]))
            if linker_tail == linker_head and len(linker_ids_sorted) > 1:
                linker_tail = sorted(
                    (idx for idx in linker_ids_sorted if idx != linker_head),
                    key=lambda idx: float(linker_positions[idx][2]),
                )[0]

            pull_anchor_atoms = {
                1: selected_atoms,
                2: [linker_head],
                3: [linker_tail],
            }
            linker_pairs.append(
                {
                    "selected_atoms": selected_atoms,
                    "instance_atoms": linker_ids_sorted,
                    "linker_head": linker_head,
                    "linker_tail": linker_tail,
                }
            )
            linker_pull_enabled = True
            print(f"  → Anchor_2 (linker head): atom {linker_head}")
            print(f"  → Anchor_3 (linker tail): atom {linker_tail}")

    if ads_mode:
        if anchor_atoms:
            print("ℹ Adsorption mode enabled: anchor restraints/pulls will be omitted.")
        pull_anchor_atoms = {}
        linker_pull_enabled = False

    # ===============================================================
    # Create folder structure
    # ===============================================================
    topo_dir = output_root / "0_topology"
    mdp_dir = output_root / "1_mdp"
    sys_dir = output_root / "2_system"

    ensure_dir(topo_dir)
    ensure_dir(mdp_dir)
    ensure_dir(sys_dir)

    shutil.copy(input_gro, sys_dir / "system.gro")

    # ===============================================================
    # Copy Martini files (UNCHANGED)
    # ===============================================================
    pkg_dir = Path(martinisurf.__file__).parent
    water_template_name = (
        POLARIZABLE_WATER_TEMPLATE if args.polarizable_water else STANDARD_WATER_TEMPLATE
    )
    water_template = pkg_dir / "system_templates" / water_template_name
    if water_template.exists():
        shutil.copy(water_template, sys_dir / water_template_name)
        print(f"✔ Copied {water_template_name} template into 2_system")
    elif args.polarizable_water:
        raise FileNotFoundError(
            f"❌ Polarizable-water template not found at {water_template}."
        )
    else:
        print(f"ℹ {water_template_name} template not found at {water_template}. Skipping solvent template copy.")

    src_itp_dir = pkg_dir / "system_itp"
    dst_itp_dir = topo_dir / "system_itp"
    ensure_dir(dst_itp_dir)

    legacy_itp_dir = output_root / "system_itp"
    if legacy_itp_dir.exists():
        for p in legacy_itp_dir.glob("*.itp"):
            dst = dst_itp_dir / p.name
            if not dst.exists():
                shutil.copy(p, dst)

    if is_dna:
        dna_ff = _select_dna_forcefield_name(
            src_itp_dir,
            polarizable_water=args.polarizable_water,
        )
        required = ["martini_v2.0_ions.itp"]
        if dna_ff:
            required.insert(0, dna_ff)
    else:
        required = [
            "martini_v3.0.0.itp",
            "martini_v3.0.0_solvents_v1.itp",
            "martini_v3.0.0_ions_v1.itp",
        ]

    for fname in required:
        src = src_itp_dir / fname
        if src.exists():
            shutil.copy(src, dst_itp_dir / fname)
    ff_moltypes: set[str] = set()
    for fname in required:
        ff_moltypes.update(_read_itp_moleculetype_set(dst_itp_dir / fname))

    if linker_mode and not (dst_itp_dir / args.linker_itp_name).exists():
        raise FileNotFoundError(
            "❌ Linker mode requested but "
            f"0_topology/system_itp/{args.linker_itp_name} is missing."
        )
    if args.cofactor_count > 0 and args.cofactor_itp_name and not (dst_itp_dir / args.cofactor_itp_name).exists():
        raise FileNotFoundError(
            "❌ Cofactor mode requested but "
            f"0_topology/system_itp/{args.cofactor_itp_name} is missing."
        )
    if args.substrate_count > 0 and args.substrate_itp_name and not (dst_itp_dir / args.substrate_itp_name).exists():
        raise FileNotFoundError(
            "❌ Substrate mode requested but "
            f"0_topology/system_itp/{args.substrate_itp_name} is missing."
        )

    # Normalize martinize POSRES_FC macros into numeric constants for robust grompp parsing.
    posres_fc_replaced = 0
    for itp_file in sorted(dst_itp_dir.glob("*.itp")):
        if _materialize_posres_fc_in_itp(itp_file):
            posres_fc_replaced += 1
    if posres_fc_replaced > 0:
        print(f"✔ Normalized POSRES_FC macros in {posres_fc_replaced} ITP file(s)")

    linker_total = 0
    linker_moltype_name: str | None = None
    if linker_mode and args.linker_resname:
        linker_atom_count = sum(
            1 for a in u.atoms if str(a.resname).strip() == args.linker_resname
        )
        if args.linker_size and args.linker_size > 0:
            linker_total = linker_atom_count // args.linker_size
        elif linker_atom_count > 0:
            linker_resids = {
                int(a.resid)
                for a in u.atoms
                if str(a.resname).strip() == args.linker_resname
            }
            linker_total = len(linker_resids)

        linker_itp_path = dst_itp_dir / args.linker_itp_name
        linker_moltype_name = _read_itp_moleculetype(linker_itp_path) or linker_itp_path.stem
    cofactor_moltype_name: str | None = None
    if args.cofactor_count > 0 and args.cofactor_itp_name:
        cofactor_itp_path = dst_itp_dir / args.cofactor_itp_name
        cofactor_moltype_name = _read_itp_moleculetype(cofactor_itp_path) or cofactor_itp_path.stem
    substrate_moltype_name: str | None = args.substrate_moltype
    if args.substrate_count > 0 and args.substrate_itp_name:
        substrate_itp_path = dst_itp_dir / args.substrate_itp_name
        substrate_moltype_name = _read_itp_moleculetype(substrate_itp_path) or substrate_itp_path.stem

    include_cofactor_itp = bool(args.cofactor_count > 0 and args.cofactor_itp_name)
    include_substrate_itp = bool(args.substrate_count > 0 and args.substrate_itp_name)
    if include_cofactor_itp and cofactor_moltype_name in ff_moltypes:
        include_cofactor_itp = False
        print(f"ℹ {cofactor_moltype_name} detected in Martini FF, skipping {args.cofactor_itp_name} include.")
    if include_substrate_itp and substrate_moltype_name in ff_moltypes:
        include_substrate_itp = False
        print(f"ℹ {substrate_moltype_name} detected in Martini FF, skipping {args.substrate_itp_name} include.")
    if args.substrate_count > 0 and not include_substrate_itp:
        if not substrate_moltype_name:
            raise ValueError("❌ Could not determine substrate moleculetype.")
        if substrate_moltype_name not in ff_moltypes:
            raise FileNotFoundError(
                "❌ Substrate topology is not available locally and "
                f"{substrate_moltype_name} was not found in the Martini force-field includes."
            )

    # ===============================================================
    # Detect molecule ITP (UNCHANGED)
    # ===============================================================
    if is_dna:
        possible_itps = [
            p for p in dst_itp_dir.glob("*.itp")
            if not p.name.startswith("martini_")
            and p.name != "surface.itp"
            and p.name != args.linker_itp_name
            and p.name != (args.cofactor_itp_name or "")
            and p.name != (args.substrate_itp_name or "")
        ]
        if not possible_itps:
            raise FileNotFoundError("❌ No DNA molecule ITP found in system_itp.")
        mol_itp = possible_itps[0]
        moltype = _read_itp_moleculetype(mol_itp) or mol_itp.stem
    else:
        if not args.moltype:
            raise ValueError("❌ --moltype required for protein")
        moltype = args.moltype
        mol_itp = dst_itp_dir / f"{moltype}.itp"
        if not mol_itp.exists():
            fallback = dst_itp_dir / "Active.itp"
            if fallback.exists():
                shutil.copy(fallback, mol_itp)
            else:
                candidates = [
                    p for p in dst_itp_dir.glob("*.itp")
                    if not p.name.startswith("martini_")
                    and p.name not in {
                        "surface.itp",
                        "go_atomtypes.itp",
                        "go_nbparams.itp",
                        args.linker_itp_name,
                        (args.cofactor_itp_name or ""),
                        (args.substrate_itp_name or ""),
                    }
                    and not p.name.endswith("_anchor.itp")
                ]
                if candidates:
                    src = sorted(candidates, key=lambda p: p.stat().st_size, reverse=True)[0]
                    _write_itp_with_moleculetype(src, mol_itp, moltype)
                else:
                    mol_itp.write_text(
                        "[ moleculetype ]\n"
                        f"{moltype} 1\n\n"
                        "[ atoms ]\n"
                        "1  C1  1  MOL  C1  1  0.0\n"
                    )

    # ===============================================================
    # Restrained anchor ITP (skip in adsorption mode)
    # ===============================================================
    active_anchor = dst_itp_dir / f"{mol_itp.stem}_anchor.itp"
    if not ads_mode:
        new_lines = []
        inside_posres = False

        with open(mol_itp) as fin:
            for line in fin:
                if "[ position_restraints" in line:
                    inside_posres = True
                    continue
                if inside_posres:
                    if line.strip().startswith("["):
                        inside_posres = False
                        new_lines.append(line)
                    continue
                new_lines.append(line)

        new_lines.append("\n[ position_restraints ]\n#ifdef POSRES\n")

        for _, atoms in anchor_atoms.items():
            for atom in atoms:
                new_lines.append(f"{atom} 1 1000 1000 0\n")

        new_lines.append("#endif\n")
        active_anchor.write_text("".join(new_lines))

    # ===============================================================
    # Copy MDP templates
    # ===============================================================
    print("• Copying MDP templates ...")

    mdp_pkg = pkg_dir / "mdp_templates"

    if is_dna:
        mdp_files = {
            "minimization_dna.mdp": "minimization_dna.mdp",
            "nvt_dna.mdp": "nvt_dna.mdp",
            "npt_dna.mdp": "npt_dna.mdp",
            "deposition_dna.mdp": "deposition_dna.mdp",
            "production_dna.mdp": "production_dna.mdp",
        }
    else:
        mdp_files = {
            "minimization.mdp": "minimization.mdp",
            "nvt.mdp": "nvt.mdp",
            "npt.mdp": "npt.mdp",
            "deposition.mdp": "deposition.mdp",
            "production.mdp": "production.mdp",
        }

    if ads_mode:
        allowed_ads = {
            "minimization.mdp",
            "nvt.mdp",
            "npt.mdp",
            "production.mdp",
            "minimization_dna.mdp",
            "nvt_dna.mdp",
            "npt_dna.mdp",
            "production_dna.mdp",
        }
        mdp_files = {src_name: dst_name for src_name, dst_name in mdp_files.items() if dst_name in allowed_ads}

    surface_itp_path = dst_itp_dir / "surface.itp"
    surface_moltype = _read_itp_moleculetype(surface_itp_path) or "SRF"
    surface_resname = _read_itp_first_atoms_resname(surface_itp_path) or surface_moltype

    # copy MDPs
    for src_name, dst_name in mdp_files.items():
        src = mdp_pkg / src_name
        if src.exists():
            rewrite_pull = dst_name in {
                "nvt.mdp",
                "npt.mdp",
                "nvt_dna.mdp",
                "npt_dna.mdp",
                "deposition.mdp",
                "production.mdp",
                "deposition_dna.mdp",
                "production_dna.mdp",
            }
            include_pull_init = rewrite_pull
            dst_path = mdp_dir / dst_name
            write_custom_mdp(
                src=src,
                dst=dst_path,
                anchor_count=len(pull_anchor_atoms),
                is_dna=is_dna,
                surface_moltype=surface_moltype,
                linker_pull=linker_pull_enabled,
                linker_count=(len(pull_anchor_atoms) // 3) if linker_pull_enabled else 1,
                linker_pull_init_prot_nm=args.linker_pull_init_prot,
                linker_pull_init_surf_nm=args.linker_pull_init_surf,
                rewrite_pull=rewrite_pull,
                include_pull_init=include_pull_init,
                dna_linker_bonded=bool(is_dna and linker_pull_enabled),
            )
            if is_dna and args.polarizable_water:
                _rewrite_mdp_for_polarizable_water(dst_path)
            print(f"  ✔ {src_name} → {dst_name}")
        else:
            print(f"⚠ Missing MDP template: {src_name}")

    anchor_moltype = mol_itp.stem if is_dna else moltype
    restrained_mol_itp_name = f"{anchor_moltype}_anchor.itp"
    if ads_mode:
        restrained_mol_itp_name = mol_itp.name
    restrained_linker_itp_name: str | None = None
    if linker_mode:
        # In linker mode, do not restrain DNA/Protein by default.
        restrained_mol_itp_name = mol_itp.name
        linker_anchor_atoms = _infer_linker_restrained_atom_ids(
            pull_anchor_atoms=pull_anchor_atoms,
            linker_atom_ids=linker_atom_ids,
            linker_size=args.linker_size,
        )
        linker_itp_path = dst_itp_dir / args.linker_itp_name
        if linker_anchor_atoms and linker_itp_path.exists():
            restrained_linker_itp_name = f"{Path(args.linker_itp_name).stem}_anchor.itp"
            _rewrite_itp_with_posres(
                src_itp=linker_itp_path,
                dst_itp=dst_itp_dir / restrained_linker_itp_name,
                posres_atom_ids=linker_anchor_atoms,
            )
            print(
                "✔ Linker restrained topology generated: "
                f"{restrained_linker_itp_name} (atoms {linker_anchor_atoms})"
            )
        else:
            print("ℹ Linker restrained topology could not be generated; skipping linker POSRES include.")

    surface_itp_path = dst_itp_dir / "surface.itp"
    surface_moltype = _read_itp_moleculetype(surface_itp_path) or "SRF"
    surface_itp_atoms = _count_itp_atoms(surface_itp_path)
    surface_atom_total = sum(
        1 for a in u.atoms if str(a.resname).strip() == surface_moltype
    )
    surface_count = 1
    if surface_itp_atoms > 0 and surface_atom_total > 0:
        surface_count = max(1, surface_atom_total // surface_itp_atoms)
        if surface_atom_total % surface_itp_atoms != 0:
            print(
                "⚠ Surface atom count is not an exact multiple of surface.itp atoms-per-molecule. "
                f"Using floor count: {surface_count}."
            )
    elif surface_atom_total > 0:
        surface_count = surface_atom_total

    intermolecular_itp_name: str | None = None
    top_mol_itp_name = mol_itp.name
    top_anchor_itp_name = restrained_mol_itp_name
    top_use_linker = linker_mode
    top_linker_count = linker_total
    top_restrained_linker_itp_name = restrained_linker_itp_name

    if is_dna and linker_mode and linker_pairs:
        bond_length_nm = args.linker_pull_init_prot if args.linker_pull_init_prot else 0.47
        merged_name = f"{mol_itp.stem}_linker.itp"
        n_bonds, n_angles = _merge_dna_linker_itp(
            dst_itp_path=dst_itp_dir / merged_name,
            dna_itp_path=mol_itp,
            linker_itp_path=dst_itp_dir / args.linker_itp_name,
            universe=u,
            linker_pairs=linker_pairs,
            linker_resname=args.linker_resname or "LNK",
            bond_length_nm=bond_length_nm,
            bond_force=1250.0,
            angle_deg=180.0,
            angle_force=20.0,
        )
        if n_bonds > 0:
            # DNA+linker becomes one molecule topology; linker no longer appears as separate molecule.
            top_mol_itp_name = merged_name
            top_anchor_itp_name = merged_name
            top_use_linker = False
            top_linker_count = 0
            top_restrained_linker_itp_name = None
            print(
                "✔ DNA-linker bonded coupling generated in merged topology: "
                f"{merged_name} ({n_bonds} bond(s), {n_angles} angle(s))"
            )
        else:
            print(
                "⚠ DNA-linker bonded coupling skipped: no BB1/BB2/BB3 bead found "
                "in selected anchor residues."
            )

    write_top_files(
        topo_dir=topo_dir,
        dst_itp_dir=dst_itp_dir,
        moltype=moltype,
        mol_itp_name=top_mol_itp_name,
        anchor_itp_name=top_anchor_itp_name,
        is_dna=is_dna,
        use_linker=top_use_linker,
        go_model=args.go_model,
        polarizable_water=args.polarizable_water,
        linker_itp_name=args.linker_itp_name,
        restrained_linker_itp_name=top_restrained_linker_itp_name,
        linker_moltype=linker_moltype_name,
        linker_count=top_linker_count,
        cofactor_itp_name=args.cofactor_itp_name if include_cofactor_itp else None,
        cofactor_moltype=cofactor_moltype_name,
        cofactor_count=args.cofactor_count,
        substrate_itp_name=args.substrate_itp_name if include_substrate_itp else None,
        substrate_moltype=substrate_moltype_name,
        substrate_count=args.substrate_count,
        surface_moltype=surface_moltype,
        surface_count=surface_count,
        intermolecular_itp_name=intermolecular_itp_name,
    )

    # ===============================================================
    # INDEX
    # ===============================================================
    groups_for_index = pull_anchor_atoms if pull_anchor_atoms else anchor_atoms
    if ads_mode:
        groups_for_index = {}
    custom_groups: Dict[str, List[int]] = {}

    # Lowercase alias used by current mdp templates (tc-grps = system).
    custom_groups["system"] = [int(a.index + 1) for a in u.atoms]

    def _collect_resname_group(*resnames: str) -> List[int]:
        wanted = {str(r).strip() for r in resnames if str(r).strip()}
        return [int(a.index + 1) for a in u.atoms if str(a.resname).strip() in wanted]

    def _iter_extra_resnames() -> List[str]:
        protein_resnames = {
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
            "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
            "THR", "TRP", "TYR", "VAL",
        }
        ignored = {
            "",
            surface_moltype,
            surface_resname,
            "W",
            "WF",
            "PW",
            "SOL",
            "NA",
            "CL",
            "K",
            "CA",
            "MG",
            "ZN",
            "LI",
            "RB",
            "CS",
            "BA",
            "SR",
            "F",
            "BR",
            "I",
            "DA",
            "DC",
            "DG",
            "DT",
        } | protein_resnames
        seen: set[str] = set()
        ordered: List[str] = []
        for atom in u.atoms:
            resname = str(atom.resname).strip()
            if resname in ignored or resname in seen:
                continue
            seen.add(resname)
            ordered.append(resname)
        return ordered

    # Surface freeze/pull groups referenced by mdp templates.
    surface_atoms = _collect_resname_group(surface_moltype, surface_resname)
    if surface_atoms:
        custom_groups[surface_moltype] = surface_atoms
        if surface_moltype != "SRF":
            custom_groups["SRF"] = surface_atoms

    # Common optional groups used in workflows and useful for downstream grompp calls.
    for group_name, resnames in (
        ("W", ("W",)),
        ("WF", ("WF",)),
        ("IONS", ("NA", "CL", "K", "CA", "MG", "ZN", "LI", "RB", "CS", "BA", "SR", "F", "BR", "I")),
        ("LINKER", (linker_moltype_name,) if linker_mode and linker_moltype_name else ()),
    ):
        atoms = _collect_resname_group(*resnames)
        if atoms:
            custom_groups[group_name] = atoms

    if is_dna:
        dna_atoms = _collect_resname_group("DA", "DC", "DG", "DT")
        if dna_atoms:
            custom_groups["DNA"] = dna_atoms
    else:
        protein_resnames = (
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
            "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
            "THR", "TRP", "TYR", "VAL",
        )
        protein_atoms = _collect_resname_group(*protein_resnames)
        if protein_atoms:
            custom_groups["Protein"] = protein_atoms

    for resname in _iter_extra_resnames():
        atoms = _collect_resname_group(resname)
        if atoms:
            custom_groups[resname] = atoms

    if groups_for_index or custom_groups:
        with open(topo_dir / "index.ndx", "w") as ndx:
            for group_name, atoms in custom_groups.items():
                if atoms:
                    ndx.write(f"\n[ {group_name} ]\n")
                    write_list(atoms, ndx)
            for gid, atoms in groups_for_index.items():
                ndx.write(f"\n[ Anchor_{gid} ]\n")
                write_list(atoms, ndx)

    if args.go_model:
        print("✔ GoMartini system built")
    else:
        print("✔ MartiniSurf system built")


if __name__ == "__main__":
    main()
