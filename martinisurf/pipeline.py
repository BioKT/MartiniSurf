#!/usr/bin/env python3
"""
MartiniSurf Full Pipeline (Protein + DNA Compatible)

Stable architecture:
• Protein → martinize2
• DNA     → martinize-dna.py
• Classical anchor mode
• Multi-linker mode
• Optional random surface linkers
"""

import argparse
import importlib.util
import os
import random
import re
import shutil
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Any, Optional

from martinisurf.utils.pdb_generation import load_clean_pdb
from martinisurf.utils.pdb_to_gro import pdb_to_gro


def _read_gro_first_resname(gro_path: str) -> str | None:
    with open(gro_path, "r") as fh:
        lines = fh.readlines()
    for line in lines[2:-1]:
        if len(line) >= 10:
            return line[5:10].strip() or None
    return None


def _read_gro_first_atomname(gro_path: str) -> str | None:
    with open(gro_path, "r") as fh:
        lines = fh.readlines()
    for line in lines[2:-1]:
        if len(line) >= 15:
            return line[10:15].strip() or None
    return None


def _read_gro_last_atomname(gro_path: str) -> str | None:
    with open(gro_path, "r") as fh:
        lines = fh.readlines()[2:-1]
    for line in reversed(lines):
        if len(line) >= 15:
            name = line[10:15].strip()
            if name:
                return name
    return None


def _read_gro_atomname_for_resid(gro_path: str, resid: int) -> str | None:
    with open(gro_path, "r") as fh:
        lines = fh.readlines()[2:-1]
    for line in lines:
        if len(line) < 15:
            continue
        try:
            gro_resid = int(line[0:5])
        except ValueError:
            continue
        if gro_resid == resid:
            name = line[10:15].strip()
            if name:
                return name
    return None


def _read_gro_atom_count(gro_path: str) -> int | None:
    with open(gro_path, "r") as fh:
        lines = fh.readlines()
    if len(lines) < 2:
        return None
    try:
        return int(lines[1].strip())
    except ValueError:
        return None


def _write_minimal_surface_itp(itp_path: Path, resname: str, bead: str, charge: float = 0.0) -> None:
    clean_res = (resname or "SRF").strip()[:6] or "SRF"
    clean_bead = (bead or "C1").strip()[:6] or "C1"
    atom_name = clean_bead
    itp_path.write_text(
        ";;;;;; Minimal surface topology (auto-generated fallback)\n\n"
        "[ moleculetype ]\n"
        "; molname nrexcl\n"
        f"  {clean_res}        1\n\n"
        "[ atoms ]\n"
        f"  1   {clean_bead:<6}   1   {clean_res:<4}   {atom_name:<5} 1     {float(charge):.3f}\n"
    )


def _read_gro_records(gro_path: str) -> tuple[str, list[dict], list[float]]:
    with open(gro_path, "r") as fh:
        lines = fh.readlines()
    if len(lines) < 3:
        raise ValueError(f"Invalid GRO file: {gro_path}")

    title = lines[0].rstrip("\n")
    atom_lines = lines[2:-1]
    box_line = lines[-1].strip()
    box_tokens = box_line.split()
    if len(box_tokens) < 3:
        raise ValueError(f"Invalid GRO box line in: {gro_path}")
    box = [float(box_tokens[0]), float(box_tokens[1]), float(box_tokens[2])]

    records: list[dict] = []
    for line in atom_lines:
        if len(line) < 44:
            continue
        try:
            records.append(
                {
                    "resid": int(line[0:5]),
                    "resname": line[5:10].strip() or "SUB",
                    "atomname": line[10:15].strip() or "C1",
                    "atomid": int(line[15:20]),
                    "x": float(line[20:28]),
                    "y": float(line[28:36]),
                    "z": float(line[36:44]),
                }
            )
        except ValueError:
            continue
    return title, records, box


def _write_gro_records(gro_path: str, title: str, records: list[dict], box: list[float]) -> None:
    with open(gro_path, "w") as fh:
        fh.write(f"{title}\n")
        fh.write(f"{len(records):5d}\n")
        for rec in records:
            resid = int(rec["resid"]) % 100000
            atomid = int(rec["atomid"]) % 100000
            fh.write(
                f"{resid:5d}"
                f"{str(rec['resname'])[:5]:<5}"
                f"{str(rec['atomname'])[:5]:>5}"
                f"{atomid:5d}"
                f"{float(rec['x']):8.3f}{float(rec['y']):8.3f}{float(rec['z']):8.3f}\n"
            )
        fh.write(f"{box[0]:10.5f}{box[1]:10.5f}{box[2]:10.5f}\n")


def _resolve_sidecar_itp(gro_path: str, explicit_itp: str | None, label: str) -> Path:
    if explicit_itp:
        itp = Path(explicit_itp)
    else:
        itp = Path(gro_path).with_suffix(".itp")
    if not itp.exists():
        raise FileNotFoundError(
            f"{label} ITP not found: {itp}. "
            f"Provide --{label.lower()}-itp or place it next to the GRO file."
        )
    return itp


def _append_random_substrates_to_gro(
    system_gro: str,
    substrate_gro: str,
    substrate_count: int,
    min_distance_nm: float = 0.20,
    max_attempts_per_copy: int = 500,
) -> None:
    if substrate_count <= 0:
        return

    title, system_records, box = _read_gro_records(system_gro)
    _, substrate_records, _ = _read_gro_records(substrate_gro)
    if not substrate_records:
        raise ValueError(f"Substrate GRO has no atoms: {substrate_gro}")

    sx = [r["x"] for r in substrate_records]
    sy = [r["y"] for r in substrate_records]
    sz = [r["z"] for r in substrate_records]
    cx = sum(sx) / len(sx)
    cy = sum(sy) / len(sy)
    cz = sum(sz) / len(sz)
    rel = [(r["x"] - cx, r["y"] - cy, r["z"] - cz) for r in substrate_records]

    min_rx = min(v[0] for v in rel)
    max_rx = max(v[0] for v in rel)
    min_ry = min(v[1] for v in rel)
    max_ry = max(v[1] for v in rel)
    min_rz = min(v[2] for v in rel)
    max_rz = max(v[2] for v in rel)

    margin = 0.05
    x_low, x_high = -min_rx + margin, box[0] - max_rx - margin
    y_low, y_high = -min_ry + margin, box[1] - max_ry - margin
    z_low, z_high = -min_rz + margin, box[2] - max_rz - margin
    if x_low > x_high or y_low > y_high or z_low > z_high:
        raise ValueError(
            "Substrate does not fit inside the simulation box. "
            "Use a larger box or a smaller substrate."
        )

    existing_xyz = [(r["x"], r["y"], r["z"]) for r in system_records]
    min_d2 = min_distance_nm * min_distance_nm
    next_atomid = max((int(r["atomid"]) for r in system_records), default=0) + 1
    next_resid = max((int(r["resid"]) for r in system_records), default=0) + 1

    original_resids: list[int] = []
    for rec in substrate_records:
        resid = int(rec["resid"])
        if resid not in original_resids:
            original_resids.append(resid)

    rng = random.Random()

    for copy_idx in range(substrate_count):
        placed = False
        for _ in range(max_attempts_per_copy):
            tx = rng.uniform(x_low, x_high)
            ty = rng.uniform(y_low, y_high)
            tz = rng.uniform(z_low, z_high)
            candidate_xyz = [(dx + tx, dy + ty, dz + tz) for dx, dy, dz in rel]

            too_close = False
            for px, py, pz in candidate_xyz:
                for ex, ey, ez in existing_xyz:
                    dx = px - ex
                    dy = py - ey
                    dz = pz - ez
                    if (dx * dx + dy * dy + dz * dz) < min_d2:
                        too_close = True
                        break
                if too_close:
                    break
            if too_close:
                continue

            resid_map = {rid: next_resid + i for i, rid in enumerate(original_resids)}
            next_resid += len(original_resids)

            for template, (x, y, z) in zip(substrate_records, candidate_xyz):
                system_records.append(
                    {
                        "resid": resid_map[int(template["resid"])],
                        "resname": template["resname"],
                        "atomname": template["atomname"],
                        "atomid": next_atomid,
                        "x": x,
                        "y": y,
                        "z": z,
                    }
                )
                next_atomid += 1
            existing_xyz.extend(candidate_xyz)
            placed = True
            break

        if not placed:
            raise RuntimeError(
                f"Could not place substrate copy {copy_idx + 1}/{substrate_count} "
                "without overlaps. Try fewer substrates or a larger box."
            )

    _write_gro_records(system_gro, title, system_records, box)


def _bead_size_class(bead_name: Optional[str]) -> str:
    if not bead_name:
        return "R"
    name = bead_name.strip().upper()
    if name.startswith("T"):
        return "T"
    if name.startswith("S"):
        return "S"
    return "R"


def _sigma_nm(is_dna: bool, ff_name: str, class_a: str, class_b: str) -> float:
    if is_dna or "martini2" in ff_name.lower():
        return 0.47

    pair = tuple(sorted((class_a, class_b)))
    table = {
        ("R", "R"): 0.47,
        ("R", "S"): 0.43,
        ("R", "T"): 0.395,
        ("S", "S"): 0.41,
        ("S", "T"): 0.365,
        ("T", "T"): 0.34,
    }
    return table.get(pair, 0.47)


def _normalize_merge_groups(merge_values: Optional[list[str]]) -> list[str]:
    if not merge_values:
        return []
    groups: list[str] = []
    for raw in merge_values:
        if raw is None:
            continue
        item = str(raw).strip()
        if not item:
            continue
        groups.append(item)
    return groups


def _parse_simple_yaml_value(raw: str) -> Any:
    text = raw.strip()
    if not text:
        return ""
    if (text.startswith('"') and text.endswith('"')) or (text.startswith("'") and text.endswith("'")):
        return text[1:-1]
    if text.lower() in {"true", "false"}:
        return text.lower() == "true"
    if text.startswith("[") and text.endswith("]"):
        inner = text[1:-1].strip()
        if not inner:
            return []
        return [_parse_simple_yaml_value(item.strip()) for item in inner.split(",")]
    if re.fullmatch(r"[+-]?\d+", text):
        return int(text)
    if re.fullmatch(r"[+-]?\d+\.\d*", text):
        return float(text)
    return text


def _load_simple_yaml(path: Path) -> dict[str, Any]:
    data: dict[str, Any] = {}
    current_section: str | None = None

    for raw in path.read_text().splitlines():
        line = raw.split("#", 1)[0].rstrip()
        if not line.strip():
            continue
        indent = len(line) - len(line.lstrip(" "))
        stripped = line.strip()
        if ":" not in stripped:
            raise ValueError(f"Invalid YAML line in {path}: {raw}")
        key, value = stripped.split(":", 1)
        key = key.strip()
        value = value.strip()

        if indent == 0:
            if value == "":
                data[key] = {}
                current_section = key
            else:
                data[key] = _parse_simple_yaml_value(value)
                current_section = None
        else:
            if not current_section:
                raise ValueError(f"Invalid nested YAML line in {path}: {raw}")
            section = data.get(current_section)
            if not isinstance(section, dict):
                raise ValueError(f"Invalid YAML section in {path}: {current_section}")
            section[key] = _parse_simple_yaml_value(value)

    return data


def _read_itp_moleculetype_name(itp_path: Path) -> str | None:
    if not itp_path.exists():
        return None
    in_moleculetype = False
    for raw in itp_path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith("[") and "moleculetype" in line.lower():
            in_moleculetype = True
            continue
        if in_moleculetype:
            if line.startswith("["):
                break
            return line.split()[0]
    return None


def _load_pre_cg_complex_config(config_path: Path) -> dict[str, Any]:
    if not config_path.exists():
        raise FileNotFoundError(f"Complex config not found: {config_path}")

    cfg = _load_simple_yaml(config_path)
    if str(cfg.get("mode", "")).strip() != "pre_cg_complex":
        raise ValueError("complex_config.yaml must define mode: pre_cg_complex")

    input_dir = config_path.parent
    complex_gro_name = str(cfg.get("complex_gro", "")).strip()
    if not complex_gro_name:
        raise ValueError("complex_config.yaml missing required key: complex_gro")
    complex_gro = (input_dir / complex_gro_name).resolve()
    if not complex_gro.exists():
        raise FileNotFoundError(f"Complex GRO not found: {complex_gro}")

    protein = cfg.get("protein")
    if not isinstance(protein, dict):
        raise ValueError("complex_config.yaml missing required section: protein")
    protein_molname = str(protein.get("molname", "")).strip()
    if not protein_molname:
        raise ValueError("complex_config.yaml missing required key: protein.molname")
    reference_pdb_name = str(protein.get("reference_pdb", "")).strip()
    reference_pdb: Path | None = None
    if reference_pdb_name:
        reference_pdb = (input_dir / reference_pdb_name).resolve()
        if not reference_pdb.exists():
            raise FileNotFoundError(f"Protein reference PDB not found: {reference_pdb}")
    anchor_groups_cfg = protein.get("anchor_groups")
    anchor_groups: list[list[int]] = []
    anchor_landmark_mode = "residue"
    if isinstance(anchor_groups_cfg, list) and anchor_groups_cfg:
        anchor_landmark_mode = "group"
        raw_anchor_groups: list[list[str]] = []
        has_chain_syntax = False
        for raw_group in anchor_groups_cfg:
            if not isinstance(raw_group, str):
                raise ValueError("protein.anchor_groups entries must be strings like '1 8 10 11'")
            parts = raw_group.split()
            if len(parts) < 2:
                raise ValueError(
                    "Each protein.anchor_groups entry must include group id and at least one residue "
                    "(example: '1 8 10 11')"
                )
            raw_anchor_groups.append(parts)
            try:
                int(parts[0])
            except (TypeError, ValueError):
                has_chain_syntax = True

        if has_chain_syntax:
            if reference_pdb is None:
                raise ValueError(
                    "protein.anchor_groups uses chain-based syntax, so protein.reference_pdb "
                    "must point to the source PDB used to build the pre-CG complex."
                )
            anchor_groups = _normalize_cli_residue_groups(
                raw_anchor_groups,
                reference_pdb,
                "protein.anchor_groups",
            )
        else:
            for parts in raw_anchor_groups:
                try:
                    parsed = [int(x) for x in parts]
                except (TypeError, ValueError) as exc:
                    raise ValueError(
                        "protein.anchor_groups entries must contain integers only "
                        "(example: '2 1025 1027 1028')"
                    ) from exc
                anchor_groups.append(parsed)
    elif "orient_by_residues" in protein:
        orient_by = protein.get("orient_by_residues")
        if not isinstance(orient_by, list) or not orient_by:
            raise ValueError(
                "complex_config.yaml requires either protein.anchor_groups "
                "or protein.orient_by_residues as a non-empty list"
            )
        try:
            orient_residues = [int(x) for x in orient_by]
        except (TypeError, ValueError) as exc:
            raise ValueError("protein.orient_by_residues must be a list of integers") from exc
        anchor_groups = [[1] + orient_residues]
    else:
        anchor_groups = []
    # pre_cg_complex defaults: use low-Z balancing unless explicitly disabled by the user.
    balance_low_z_raw = protein.get("balance_low_z", True)
    if not isinstance(balance_low_z_raw, bool):
        raise ValueError("protein.balance_low_z must be true or false")
    balance_low_z = bool(balance_low_z_raw)
    balance_low_z_fraction_raw = protein.get("balance_low_z_fraction", 0.2)
    if not isinstance(balance_low_z_fraction_raw, (int, float)):
        raise ValueError("protein.balance_low_z_fraction must be a number in (0, 1]")
    balance_low_z_fraction = float(balance_low_z_fraction_raw)
    if not (0.0 < balance_low_z_fraction <= 1.0):
        raise ValueError("protein.balance_low_z_fraction must be in the interval (0, 1]")

    cofactor = cfg.get("cofactor")
    if not isinstance(cofactor, dict):
        raise ValueError("complex_config.yaml missing required section: cofactor")
    cofactor_molname = str(cofactor.get("molname", "")).strip()
    cofactor_itp_name = str(cofactor.get("itp", "")).strip()
    cofactor_count = cofactor.get("count")
    if not cofactor_molname:
        raise ValueError("complex_config.yaml missing required key: cofactor.molname")
    if not cofactor_itp_name:
        raise ValueError("complex_config.yaml missing required key: cofactor.itp")
    if not isinstance(cofactor_count, int) or cofactor_count < 1:
        raise ValueError("complex_config.yaml requires cofactor.count as an integer >= 1")

    topology = cfg.get("topology")
    if not isinstance(topology, dict):
        raise ValueError("complex_config.yaml missing required section: topology")
    protein_itp_name = str(topology.get("protein_itp", "")).strip()
    include_go = bool(topology.get("include_go", False))
    go_glob = str(topology.get("go_files_glob", "go_*")).strip() or "go_*"
    if not protein_itp_name:
        raise ValueError("complex_config.yaml missing required key: topology.protein_itp")

    protein_itp = (input_dir / protein_itp_name).resolve()
    if not protein_itp.exists():
        raise FileNotFoundError(f"Protein ITP not found: {protein_itp}")
    cofactor_itp = (input_dir / cofactor_itp_name).resolve()
    if not cofactor_itp.exists():
        raise FileNotFoundError(f"Cofactor ITP not found: {cofactor_itp}")

    protein_itp_moltype = _read_itp_moleculetype_name(protein_itp)
    if protein_itp_moltype and protein_itp_moltype != protein_molname:
        raise ValueError(
            f"protein.molname ({protein_molname}) does not match moleculetype in {protein_itp.name} "
            f"({protein_itp_moltype})"
        )
    cofactor_itp_moltype = _read_itp_moleculetype_name(cofactor_itp)
    if cofactor_itp_moltype and cofactor_itp_moltype != cofactor_molname:
        raise ValueError(
            f"cofactor.molname ({cofactor_molname}) does not match moleculetype in {cofactor_itp.name} "
            f"({cofactor_itp_moltype})"
        )

    go_files = sorted(input_dir.glob(go_glob))
    if include_go and not go_files:
        raise FileNotFoundError(
            f"include_go is true but no files match '{go_glob}' in {input_dir}"
        )

    return {
        "complex_gro": complex_gro,
        "protein_molname": protein_molname,
        "anchor_groups": anchor_groups,
        "anchor_landmark_mode": anchor_landmark_mode,
        "balance_low_z": balance_low_z,
        "balance_low_z_fraction": balance_low_z_fraction,
        "reference_pdb": reference_pdb,
        "protein_itp": protein_itp,
        "cofactor_molname": cofactor_molname,
        "cofactor_itp": cofactor_itp,
        "cofactor_count": cofactor_count,
        "include_go": include_go,
        "go_files": go_files,
    }


def _build_pdb_chain_residue_map(pdb_path: Path) -> dict[tuple[str, int], int]:
    """
    Map (chain_id, local_resid) from the cleaned input PDB to the global residue ids
    later used in the CG GRO/system outputs.

    The global residue id follows the first-seen residue order in the cleaned PDB.
    This matches the numbering used by the downstream martinization/orientation flow.
    """
    residue_order: list[tuple[str, int, str]] = []
    seen_full_keys: set[tuple[str, int, str]] = set()
    ambiguous_keys: set[tuple[str, int]] = set()
    chain_resid_to_global: dict[tuple[str, int], int] = {}

    with open(pdb_path, "r") as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if len(line) < 27:
                continue

            chain_id = line[21].strip().upper()
            try:
                resid = int(line[22:26])
            except ValueError:
                continue
            insertion_code = line[26].strip().upper()

            full_key = (chain_id, resid, insertion_code)
            if full_key in seen_full_keys:
                continue

            seen_full_keys.add(full_key)
            residue_order.append(full_key)

            short_key = (chain_id, resid)
            if short_key in chain_resid_to_global:
                ambiguous_keys.add(short_key)
                continue

            chain_resid_to_global[short_key] = len(residue_order)

    if ambiguous_keys:
        preview = ", ".join(
            f"{chain or '?'}:{resid}"
            for chain, resid in sorted(ambiguous_keys)[:8]
        )
        raise ValueError(
            "Chain-based residue lookup is ambiguous because the cleaned PDB contains "
            f"multiple insertion-code variants for the same chain/residue id: {preview}"
        )

    return chain_resid_to_global


def _normalize_cli_residue_groups(
    raw_groups: list[list[str]] | None,
    pdb_path: Path,
    flag_name: str,
) -> list[list[int]]:
    """
    Accept legacy numeric groups (GROUP RESID...) or chain-based groups
    (CHAIN RESID...) and normalize them to the numeric format expected downstream.
    """
    if not raw_groups:
        return []

    chain_map = _build_pdb_chain_residue_map(pdb_path)
    normalized: list[list[int]] = []

    for idx, raw_group in enumerate(raw_groups, start=1):
        if len(raw_group) < 2:
            raise ValueError(
                f"{flag_name} requires at least two values: GROUP_OR_CHAIN RESID [RESID ...]"
            )

        head = str(raw_group[0]).strip()
        residue_tokens = raw_group[1:]
        try:
            residues = [int(token) for token in residue_tokens]
        except ValueError as exc:
            raise ValueError(
                f"{flag_name} residue ids must be integers. Received: {' '.join(map(str, raw_group))}"
            ) from exc

        try:
            group_id = int(head)
        except ValueError:
            chain_id = head.upper()
            resolved_residues: list[int] = []
            for resid in residues:
                key = (chain_id, resid)
                global_resid = chain_map.get(key)
                if global_resid is None:
                    available = sorted(
                        str(r)
                        for c, r in chain_map
                        if c == chain_id
                    )
                    preview = ", ".join(available[:12])
                    if len(available) > 12:
                        preview += ", ..."
                    raise ValueError(
                        f"{flag_name} chain-based selection could not resolve {chain_id} {resid}. "
                        f"Available residues for chain {chain_id}: [{preview}]"
                    )
                resolved_residues.append(global_resid)
            group_id = idx
            normalized.append([group_id] + resolved_residues)
        else:
            normalized.append([group_id] + residues)

    return normalized


def _validate_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    complex_mode = bool(args.complex_config)
    if not complex_mode and not args.pdb:
        parser.error("--pdb is required unless --complex-config is provided.")

    if not args.surface and (args.lx is None or args.ly is None):
        parser.error("When --surface is not provided, both --lx and --ly are required.")

    if args.linker and not args.linker_group:
        parser.error("Linker mode requires at least one --linker-group.")

    if not args.linker and not args.anchor and not complex_mode:
        parser.error("Provide --anchor in classical mode, or use --linker mode.")

    if args.anchor and args.linker:
        print("⚠ Both --anchor and --linker were provided. Linker mode will be used.")
    if args.ads_mode and args.linker:
        parser.error("--ads-mode is incompatible with linker mode.")
    if args.ads_mode and not (args.anchor or args.complex_config):
        parser.error("--ads-mode requires anchor-based orientation (--anchor or --complex-config).")

    if (args.go_eps is not None or args.go_low is not None or args.go_up is not None) and not args.go:
        print("⚠ Go parameters (--go-eps/--go-low/--go-up) were provided without --go. They will be ignored.")

    if args.substrate_count < 0:
        parser.error("--substrate-count must be >= 0.")
    if args.substrate and args.substrate_count == 0:
        args.substrate_count = 1
    if args.substrate_count > 0 and not args.substrate:
        parser.error("--substrate-count requires --substrate.")
    if args.ionize and not args.solvate:
        parser.error("--ionize requires --solvate.")
    if args.salt_conc < 0:
        parser.error("--salt-conc must be >= 0.")
    if args.solvate_radius <= 0:
        parser.error("--solvate-radius must be > 0.")
    if args.solvate_surface_clearance < 0:
        parser.error("--solvate-surface-clearance must be >= 0.")
    if args.balance_low_z_fraction is not None and not (0.0 < args.balance_low_z_fraction <= 1.0):
        parser.error("--balance-low-z-fraction must be in the interval (0, 1].")
    if args.water_gro and not Path(args.water_gro).exists():
        parser.error(f"--water-gro not found: {args.water_gro}")
    if args.freeze_water_fraction < 0 or args.freeze_water_fraction > 1:
        parser.error("--freeze-water-fraction must be in [0, 1].")
    if args.freeze_water_fraction > 0 and not args.dna:
        parser.error("--freeze-water-fraction is currently supported only with --dna.")
    if args.freeze_water_fraction > 0 and not args.solvate:
        parser.error("--freeze-water-fraction requires --solvate.")
    if complex_mode:
        if args.dna:
            parser.error("--complex-config currently supports protein systems only (no --dna).")
        if args.pdb:
            print("⚠ --pdb is ignored when --complex-config is provided.")


def _print_config_summary(args: argparse.Namespace) -> None:
    if args.complex_config:
        print("\n=== MartiniSurf Configuration ===")
        print("Mode:             Protein (pre-CG complex)")
        print(f"Complex config:   {args.complex_config}")
        orientation = "ads (from config)" if args.ads_mode else "anchor (from config)"
        print(f"Orientation:      {orientation}")
        print(f"Output:           {args.outdir}")
        print("=================================\n")
        return

    mode = "DNA" if args.dna else "Protein"
    orient_mode = "linker" if args.linker else "anchor"
    if args.ads_mode and not args.linker:
        orient_mode = "ads"
    print("\n=== MartiniSurf Configuration ===")
    print(f"Mode:            {mode}")
    print(f"PDB/Input:       {args.pdb}")
    print(f"Orientation:      {orient_mode}")
    print(f"Output:           {args.outdir}")
    if args.go and not args.dna:
        print("Go model:         enabled")
    if not args.dna:
        print(f"Max warnings:     {args.maxwarn}")
    if args.linker:
        group_count = len(args.linker_group) if args.linker_group else 0
        print(f"Linker groups:    {group_count}")
        print(f"Invert linker:    {args.invert_linker}")
    elif args.anchor:
        print(f"Anchor groups:    {len(args.anchor)}")
    if args.substrate and args.substrate_count > 0:
        print(f"Substrates:       {args.substrate_count}")
    if args.merge:
        print(f"Merge groups:     {', '.join(args.merge)}")
    print("=================================\n")


# ======================================================================
# PARSER
# ======================================================================

def build_parser():
    parser = argparse.ArgumentParser(
        prog="martinisurf",
        description=(
            "Build complete MartiniSurf systems for Protein/DNA on surfaces.\n"
            "Use ONE orientation mode: classical anchors (--anchor) or linker mode (--linker)."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=(
            "Examples:\n"
            "  Anchor mode:\n"
            "    martinisurf --pdb 1RJW --moltype Protein --lx 20 --ly 20 "
            "--anchor A 8 10 --anchor D 8 10\n"
            "  Linker mode:\n"
            "    martinisurf --dna --pdb 4C64.pdb --surface surface.gro "
            "--linker linker.gro --linker-group A 1\n"
        ),
    )

    input_group = parser.add_argument_group("Input And Molecule")
    input_group.add_argument("--pdb", help="Local PDB path, RCSB ID (4 chars), or UniProt ID (6 chars).")
    input_group.add_argument(
        "--complex-config",
        help=(
            "YAML config for pre-CG protein+cofactor workflow (skips martinization). "
            "Example: input/complex_config.yaml."
        ),
    )
    input_group.add_argument("--moltype", help="Molecule name for protein topology output.")
    input_group.add_argument(
        "--go",
        action="store_true",
        help="Enable Go model in martinize2 protein mode.",
    )
    input_group.add_argument("--ff", default="martini3001", help="Force field name for martinize2 (protein mode).")
    input_group.add_argument("--dna", action="store_true", help="Enable DNA mode (uses martinize-dna.py).")
    input_group.add_argument("--dnatype", default="ds-stiff", help="DNA type for martinize-dna.py.")
    input_group.add_argument(
        "--merge",
        action="append",
        metavar="CHAINS",
        help=(
            "Merge chains during martinization. Example: --merge A,B,C,D "
            "(can be repeated). Use --merge all to merge every chain."
        ),
    )

    martinize_group = parser.add_argument_group("Martinization Controls")
    martinize_group.add_argument("--p", choices=["none", "all", "backbone"], default="backbone", help="Position restraints selection.")
    martinize_group.add_argument("--pf", type=float, default=1000, help="Position restraints force constant.")
    martinize_group.add_argument("--maxwarn", type=int, default=0, help="Allowed martinize2 warnings before abort.")
    martinize_group.add_argument("--dssp", action="store_true", help="Use DSSP during protein martinization.")
    martinize_group.add_argument("--no-dssp", dest="dssp", action="store_false", help="Disable DSSP during protein martinization.")
    martinize_group.add_argument("--elastic", action="store_true", help="Enable elastic network.")
    martinize_group.add_argument("--ef", type=float, default=700, help="Elastic network force constant.")
    martinize_group.add_argument("--go-eps", type=float, help="Go model epsilon value for martinize2.")
    martinize_group.add_argument("--go-low", type=float, help="Go model minimum contact distance (nm) for martinize2.")
    martinize_group.add_argument("--go-up", type=float, help="Go model maximum contact distance (nm) for martinize2.")
    parser.set_defaults(dssp=True)

    surface_group = parser.add_argument_group("Surface (Required if --surface is omitted)")
    surface_group.add_argument("--surface", help="Existing surface .gro file. If omitted, a surface is generated.")
    surface_group.add_argument(
        "--surface-mode",
        choices=["2-1", "4-1"],
        default="2-1",
        help="Surface lattice mode for generated surfaces.",
    )
    surface_group.add_argument("--lx", type=float, help="Surface size in X (nm) for generated surface.")
    surface_group.add_argument("--ly", type=float, help="Surface size in Y (nm) for generated surface.")
    surface_group.add_argument("--dx", type=float, default=0.47, help="Surface bead spacing (nm).")
    surface_group.add_argument("--surface-bead", default="C1", help="Surface bead type for generated surface.")
    surface_group.add_argument("--charge", type=float, default=0.0, help="Surface bead charge for generated surface.")

    anchor_group = parser.add_argument_group("Orientation: Classical Anchor Mode")
    anchor_group.add_argument(
        "--anchor",
        nargs="+",
        action="append",
        metavar=("GROUP_OR_CHAIN", "RESID"),
        help=(
            "Anchor group: GROUP RESID [RESID ...] or CHAIN RESID [RESID ...]. "
            "Chain-based syntax is resolved from the cleaned input PDB."
        ),
    )
    anchor_group.add_argument("--dist", type=float, default=1.0, help="Anchor-to-surface target distance (nm).")
    anchor_group.add_argument(
        "--ads-mode",
        action="store_true",
        help=(
            "Adsorption mode: keeps anchor-based orientation but skips anchor/pull restraints in topology generation."
        ),
    )
    anchor_group.add_argument(
        "--balance-low-z",
        action="store_true",
        help="In two-anchor orientation, choose the roll angle that flattens the lowest-Z region.",
    )
    anchor_group.add_argument(
        "--balance-low-z-fraction",
        type=float,
        default=None,
        help="Fraction (0,1] of lowest-Z beads used by --balance-low-z (default: 0.2).",
    )

    linker_group = parser.add_argument_group("Orientation: Linker Mode")
    linker_group.add_argument("--linker", help="Linker .gro file.")
    linker_group.add_argument(
        "--linker-group",
        nargs="+",
        action="append",
        metavar=("GROUP_OR_CHAIN", "RESID"),
        help=(
            "Residue group(s) to attach linker(s): GROUP RESID [RESID ...] or "
            "CHAIN RESID [RESID ...]. Chain-based syntax is resolved from the cleaned input PDB."
        ),
    )
    linker_group.add_argument("--linker-prot-dist", type=float, help="Linker-to-protein/DNA distance (nm). Auto if omitted.")
    linker_group.add_argument("--linker-surf-dist", type=float, help="Linker-to-surface distance (nm). Auto if omitted.")
    linker_group.add_argument("--invert-linker", action="store_true", help="Reverse linker bead order before attachment.")
    linker_group.add_argument("--surface-linkers", type=int, default=0, help="Add random extra linkers on the surface.")

    substrate_group = parser.add_argument_group("Optional Random Substrate")
    substrate_group.add_argument("--substrate", help="Substrate .gro file to place randomly inside the simulation box.")
    substrate_group.add_argument("--substrate-itp", help="Substrate .itp file (optional; inferred from --substrate basename).")
    substrate_group.add_argument("--substrate-count", type=int, default=0, help="Number of substrate molecules to add randomly.")

    post_group = parser.add_argument_group("Optional Solvation And Ionization (requires GROMACS)")
    post_group.add_argument("--solvate", action="store_true", help="Run gmx solvate using MartiniSurf water.gro and produce final solvated files.")
    post_group.add_argument("--ionize", action="store_true", help="Run gmx genion after solvation and produce final ionized files.")
    post_group.add_argument("--salt-conc", type=float, default=0.15, help="Target salt concentration (M) for --ionize.")
    post_group.add_argument("--water-gro", help="Optional custom water coordinate file (.gro) for gmx solvate.")
    post_group.add_argument("--solvate-radius", type=float, default=0.21, help="Exclusion radius (nm) for gmx solvate (Martini recommended: 0.21).")
    post_group.add_argument("--solvate-surface-clearance", type=float, default=0.4, help="Remove water in a symmetric slab around the surface plane: |z - z_surface| <= clearance.")
    post_group.add_argument("--freeze-water-fraction", type=float, default=0.0, help="DNA-only: convert this fraction of W waters to WF in final outputs.")
    post_group.add_argument("--freeze-water-seed", type=int, default=42, help="Random seed used for DNA water freezing.")

    output_group = parser.add_argument_group("Output")
    output_group.add_argument("--outdir", default="Simulation_Files", help="Output directory for generated system files.")

    return parser


# ======================================================================
# RUNNER
# ======================================================================

def run(cmd, cwd=None):
    print("\n▶ Running:\n ", " ".join(cmd), "\n")
    res = subprocess.run(cmd, cwd=cwd)
    if res.returncode != 0:
        raise RuntimeError("Command failed.")
    print("✔ Done\n")


def _run_capture(cmd: list[str], cwd: Path | None = None, stdin_text: str | None = None) -> subprocess.CompletedProcess:
    return subprocess.run(
        cmd,
        cwd=cwd,
        text=True,
        input=stdin_text,
        capture_output=True,
        check=False,
    )


def _run_with_check(cmd: list[str], cwd: Path | None = None, stdin_text: str | None = None) -> subprocess.CompletedProcess:
    res = _run_capture(cmd, cwd=cwd, stdin_text=stdin_text)
    if res.returncode != 0:
        shown = " ".join(cmd)
        raise RuntimeError(
            f"Command failed ({shown})\nSTDOUT:\n{res.stdout[-4000:]}\nSTDERR:\n{res.stderr[-4000:]}"
        )
    return res


def _find_gmx_binary() -> str | None:
    for exe in ("gmx", "gmx_mpi"):
        if shutil.which(exe):
            return exe
    return None


def _write_ions_mdp(mdp_path: Path) -> None:
    mdp_path.write_text(
        "integrator = steep\n"
        "nsteps = 50\n"
        "emtol = 1000\n"
        "emstep = 0.01\n"
        "cutoff-scheme = Verlet\n"
        "nstlist = 20\n"
        "coulombtype = reaction-field\n"
        "rcoulomb = 1.1\n"
        "epsilon_r = 15\n"
        "epsilon_rf = 0\n"
        "vdwtype = cutoff\n"
        "rvdw = 1.1\n"
        "pbc = xyz\n"
    )


def _read_itp_moleculetype(itp_path: Path) -> str | None:
    if not itp_path.exists():
        return None
    in_moleculetype = False
    for raw in itp_path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith("[") and "moleculetype" in line.lower():
            in_moleculetype = True
            continue
        if in_moleculetype:
            if line.startswith("["):
                break
            return line.split()[0]
    return None


def _read_itp_atomname_for_moltype(itp_path: Path, moltype: str) -> str | None:
    if not itp_path.exists():
        return None

    current_moltype: str | None = None
    in_moleculetype = False
    in_atoms = False

    for raw in itp_path.read_text().splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line:
            continue

        if line.startswith("["):
            section = line.strip("[]").strip().lower()
            in_moleculetype = section == "moleculetype"
            in_atoms = section == "atoms"
            continue

        if in_moleculetype:
            current_moltype = line.split()[0]
            in_moleculetype = False
            continue

        if in_atoms and current_moltype == moltype:
            parts = line.split()
            if len(parts) >= 5:
                return parts[4]

    return None


def _read_itp_uniform_atomname_for_moltype(itp_path: Path, moltype: str) -> str | None:
    if not itp_path.exists():
        return None

    current_moltype: str | None = None
    in_moleculetype = False
    in_atoms = False
    atom_names: list[str] = []

    for raw in itp_path.read_text().splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line:
            continue

        if line.startswith("["):
            section = line.strip("[]").strip().lower()
            in_moleculetype = section == "moleculetype"
            in_atoms = section == "atoms"
            continue

        if in_moleculetype:
            current_moltype = line.split()[0]
            in_moleculetype = False
            atom_names = []
            continue

        if in_atoms and current_moltype == moltype:
            parts = line.split()
            if len(parts) >= 5:
                atom_names.append(parts[4])

    if not atom_names:
        return None
    unique = set(atom_names)
    if len(unique) != 1:
        return None
    return atom_names[0]


def _read_itp_moleculetype_names(itp_path: Path) -> list[str]:
    if not itp_path.exists():
        return []

    names: list[str] = []
    in_moleculetype = False
    for raw in itp_path.read_text().splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line:
            continue
        if line.startswith("["):
            section = line.strip("[]").strip().lower()
            in_moleculetype = section == "moleculetype"
            continue
        if in_moleculetype:
            names.append(line.split()[0])
            in_moleculetype = False
    return names


def _load_moltype_atoms_from_itp_dir(itp_dir: Path) -> dict[str, list[str]]:
    mol_atoms: dict[str, list[str]] = {}
    if not itp_dir.exists():
        return mol_atoms

    for itp in sorted(itp_dir.glob("*.itp")):
        current_moltype: str | None = None
        in_moleculetype = False
        in_atoms = False
        collecting_atoms: list[str] | None = None

        for raw in itp.read_text().splitlines():
            line = raw.split(";", 1)[0].strip()
            if not line:
                continue

            if line.startswith("["):
                # Finalize previously collected atoms before switching section.
                if collecting_atoms is not None and current_moltype and current_moltype not in mol_atoms:
                    mol_atoms[current_moltype] = collecting_atoms
                    collecting_atoms = None
                section = line.strip("[]").strip().lower()
                in_moleculetype = section == "moleculetype"
                in_atoms = section == "atoms"
                continue

            if in_moleculetype:
                current_moltype = line.split()[0]
                in_moleculetype = False
                collecting_atoms = []
                continue

            if in_atoms and current_moltype:
                parts = line.split()
                if len(parts) >= 5:
                    if collecting_atoms is not None:
                        collecting_atoms.append(parts[4])

        if collecting_atoms is not None and current_moltype and current_moltype not in mol_atoms:
            mol_atoms[current_moltype] = collecting_atoms

    return mol_atoms


def _validate_named_molecule_atomnames(top_dir: Path, top_path: Path, gro_path: Path) -> None:
    """
    Validate molecules whose GRO residue name equals their moleculetype
    (surface/cofactor/substrate/linker/ions/water).
    """
    if not top_path.exists() or not gro_path.exists():
        return

    molecules = _parse_molecules_entries(top_path)
    if not molecules:
        return

    mol_atoms = _load_moltype_atoms_from_itp_dir(top_dir / "system_itp")
    _, records, _ = _read_gro_records(str(gro_path))

    grouped: dict[tuple[int, str], list[dict]] = {}
    for rec in records:
        key = (int(rec["resid"]), str(rec["resname"]).strip())
        grouped.setdefault(key, []).append(rec)

    by_resname: dict[str, list[list[dict]]] = {}
    for (_, resname), atoms in grouped.items():
        by_resname.setdefault(resname, []).append(atoms)
    atomnames_by_resname: dict[str, list[str]] = {}
    for rec in records:
        res = str(rec["resname"]).strip()
        atomnames_by_resname.setdefault(res, []).append(str(rec["atomname"]).strip())

    issues: list[str] = []
    for molname, count in molecules:
        if count <= 0:
            continue
        atom_template = mol_atoms.get(molname)
        got_atomnames = atomnames_by_resname.get(molname, [])
        if not atom_template or not got_atomnames:
            continue
        expected_total = int(count) * len(atom_template)
        if len(got_atomnames) != expected_total:
            issues.append(
                f"{molname}: expected {expected_total} atoms from topology, found {len(got_atomnames)} atoms in GRO"
            )
            continue

        expected = [str(a).strip() for a in atom_template]
        if len(set(expected)) == 1:
            target = expected[0]
            for idx, got_name in enumerate(got_atomnames, start=1):
                if got_name != target:
                    issues.append(
                        f"{molname}: atom occurrence {idx} expected '{target}' got '{got_name}'"
                    )
                if len(issues) >= 200:
                    issues.append("... truncated ...")
                    break
            if len(issues) >= 200:
                break
            continue

        groups = by_resname.get(molname, [])
        if len(groups) != count:
            # For non-uniform atomname templates we need one GRO-residue per molecule
            # to validate sequence safely.
            continue

        for idx, atoms in enumerate(groups, start=1):
            got = [str(a["atomname"]).strip() for a in atoms]
            if len(got) != len(expected):
                continue
            mismatch_pos = next((i for i, (e, g) in enumerate(zip(expected, got), start=1) if e != g), None)
            if mismatch_pos is not None:
                issues.append(
                    f"{molname} molecule {idx}: atom {mismatch_pos} expected '{expected[mismatch_pos - 1]}' got '{got[mismatch_pos - 1]}'"
                )
            if len(issues) >= 200:
                issues.append("... truncated ...")
                break
        if len(issues) >= 200:
            break

    if not issues:
        return

    report_path = top_dir / "atomname_validation_report.txt"
    lines = [
        "MartiniSurf atomname validation report",
        f"Topology: {top_path}",
        f"Structure: {gro_path}",
        "",
        "Detected mismatches:",
        *issues,
    ]
    report_path.write_text("\n".join(lines) + "\n")
    raise RuntimeError(
        "Atomname mismatch detected between topology and GRO. "
        f"Review and fix names, then rerun. Report: {report_path}"
    )


def _detect_ff_ion_moltypes(top_dir: Path) -> dict[str, str]:
    itp_dir = top_dir / "system_itp"
    ion_itp_candidates = [
        itp_dir / "martini_v2.0_ions.itp",
        itp_dir / "martini_v3.0.0_ions_v1.itp",
    ]

    available: set[str] = set()
    for itp in ion_itp_candidates:
        available.update(_read_itp_moleculetype_names(itp))

    mapping = {
        "NA": "NA+"
        if "NA+" in available
        else ("NA" if "NA" in available else "NA"),
        "CL": "CL-"
        if "CL-" in available
        else ("CL" if "CL" in available else "CL"),
    }
    return mapping


def _normalize_ion_atom_names_from_itp(top_dir: Path, gro_path: Path) -> bool:
    itp_dir = top_dir / "system_itp"
    ion_itp_candidates = [
        itp_dir / "martini_v2.0_ions.itp",
        itp_dir / "martini_v3.0.0_ions_v1.itp",
    ]

    ff_ions = _detect_ff_ion_moltypes(top_dir)
    atomname_by_moltype: dict[str, str] = {}
    for itp in ion_itp_candidates:
        if not itp.exists():
            continue
        for mol in set(ff_ions.values()):
            atom_name = _read_itp_atomname_for_moltype(itp, mol)
            if atom_name:
                atomname_by_moltype[mol] = atom_name

    if not gro_path.exists():
        return False

    resname_to_target = {
        "NA": ff_ions["NA"],
        "NA+": ff_ions["NA"],
        "CL": ff_ions["CL"],
        "CL-": ff_ions["CL"],
    }

    title, records, box = _read_gro_records(str(gro_path))
    changed = False
    for rec in records:
        res = str(rec["resname"]).strip()
        target_res = resname_to_target.get(res)
        if not target_res:
            continue
        if res != target_res:
            rec["resname"] = target_res
            changed = True
        target_atom = atomname_by_moltype.get(target_res)
        if target_atom and str(rec["atomname"]).strip() != target_atom:
            rec["atomname"] = target_atom
            changed = True

    if changed:
        _write_gro_records(str(gro_path), title, records, box)
    return changed


def _normalize_uniform_atom_names_from_itp(top_dir: Path, top_path: Path, gro_path: Path) -> bool:
    if not top_path.exists() or not gro_path.exists():
        return False

    entries = _parse_molecules_entries(top_path)
    if not entries:
        return False

    water_ion_names = {
        "W", "SOL", "WF", "NA", "NA+", "CL", "CL-", "K", "CA", "MG", "ZN", "LI", "RB", "CS", "BA", "SR", "F", "BR", "I",
    }
    moltypes = [name for name, count in entries if count > 0 and name not in water_ion_names]
    if not moltypes:
        return False

    uniform_map: dict[str, str] = {}
    for itp in sorted((top_dir / "system_itp").glob("*.itp")):
        for mol in moltypes:
            if mol in uniform_map:
                continue
            atom = _read_itp_uniform_atomname_for_moltype(itp, mol)
            if atom:
                uniform_map[mol] = atom

    if not uniform_map:
        return False

    title, records, box = _read_gro_records(str(gro_path))
    changed = False
    for rec in records:
        res = str(rec["resname"]).strip()
        target_atom = uniform_map.get(res)
        if target_atom and str(rec["atomname"]).strip() != target_atom:
            rec["atomname"] = target_atom
            changed = True

    if changed:
        _write_gro_records(str(gro_path), title, records, box)
    return changed


def _count_residues_by_resname(records: list[dict], target_resnames: set[str]) -> int:
    keys = {
        (int(rec["resid"]), str(rec["resname"]).strip())
        for rec in records
        if str(rec["resname"]).strip() in target_resnames
    }
    return len(keys)


def _update_top_molecule_count(top_path: Path, molname: str, new_count: int) -> None:
    def _fmt_molecule_line(name: str, count: int) -> str:
        # Keep a stable fixed-width layout for [ molecules ] entries.
        return f"{name:<16}{int(count)}"

    lines = top_path.read_text().splitlines()
    out: list[str] = []
    in_molecules = False
    replaced = False
    for raw in lines:
        stripped = raw.strip()
        if stripped.lower() == "[ molecules ]":
            in_molecules = True
            out.append(raw)
            continue
        if in_molecules:
            if stripped.startswith("["):
                in_molecules = False
                out.append(raw)
                continue
            if stripped and not stripped.startswith(";"):
                parts = stripped.split()
                if parts and parts[0] == molname:
                    out.append(_fmt_molecule_line(molname, new_count))
                    replaced = True
                    continue
        out.append(raw)

    if not replaced:
        text = "\n".join(out).rstrip() + "\n"
        if "[ molecules ]" in text:
            text += _fmt_molecule_line(molname, new_count) + "\n"
        else:
            text += "\n[ molecules ]\n" + _fmt_molecule_line(molname, new_count) + "\n"
        top_path.write_text(text)
    else:
        top_path.write_text("\n".join(out).rstrip() + "\n")


def _remove_waters_below_surface(
    gro_path: Path,
    top_path: Path,
    surface_resname: str,
    clearance_nm: float,
    water_resnames: set[str],
    water_molname: str,
) -> None:
    title, records, box = _read_gro_records(str(gro_path))
    if not records:
        return

    surf_coords = [float(r["z"]) for r in records if str(r["resname"]).strip() == surface_resname]
    if not surf_coords:
        return
    # Use the surface plane center and remove water in a symmetric slab around it.
    z_surface = sum(surf_coords) / len(surf_coords)
    clearance = float(clearance_nm)

    # Build molecule-level groups to evaluate water COM against the slab,
    # then filter by original atom order to avoid reordering GO virtual sites.
    grouped: dict[tuple[int, str], list[dict]] = {}
    for rec in records:
        key = (int(rec["resid"]), str(rec["resname"]).strip())
        grouped.setdefault(key, []).append(rec)

    removed_keys: set[tuple[int, str]] = set()
    for key, atoms in grouped.items():
        _, resname = key
        if resname not in water_resnames:
            continue
        z_center = sum(float(a["z"]) for a in atoms) / len(atoms)
        if abs(z_center - z_surface) <= clearance:
            removed_keys.add(key)

    removed_water_res = len(removed_keys)
    kept: list[dict] = [
        rec
        for rec in records
        if (int(rec["resid"]), str(rec["resname"]).strip()) not in removed_keys
    ]

    if removed_water_res == 0:
        return

    for i, rec in enumerate(kept, start=1):
        rec["atomid"] = i
    _write_gro_records(str(gro_path), title, kept, box)

    water_res_total = _count_residues_by_resname(kept, water_resnames)
    _update_top_molecule_count(top_path, water_molname, water_res_total)
    # Silent cleanup to keep terminal output concise.


def _run_genion_with_fallback(
    gmx_bin: str,
    tpr_path: Path,
    out_gro: Path,
    top_path: Path,
    salt_conc: float,
    cwd: Path | None = None,
) -> None:
    cmd = [
        gmx_bin, "genion",
        "-s", str(tpr_path),
        "-o", str(out_gro),
        "-p", str(top_path),
        "-pname", "NA",
        "-nname", "CL",
        "-neutral",
    ]
    if salt_conc > 0:
        cmd += ["-conc", str(salt_conc)]

    candidates = ["W\n", "SOL\n"] + [f"{i}\n" for i in range(0, 50)]
    last_err = ""
    for selection in candidates:
        res = _run_capture(cmd, cwd=cwd, stdin_text=selection)
        if res.returncode == 0:
            return
        last_err = f"STDOUT:\n{res.stdout[-2000:]}\nSTDERR:\n{res.stderr[-2000:]}"

    raise RuntimeError(
        "gmx genion failed for all tested solvent selections (W/SOL/group numbers 0-49).\n"
        f"{last_err}"
    )


def _restore_non_solvent_from_reference(
    reference_gro: Path,
    target_gro: Path,
    water_resnames: set[str],
    ion_resnames: set[str] | None = None,
) -> bool:
    ions = ion_resnames or {"NA", "CL"}
    _, ref_records, _ = _read_gro_records(str(reference_gro))
    title, tgt_records, box = _read_gro_records(str(target_gro))

    def _is_solvent_or_ion(rec: dict) -> bool:
        res = str(rec["resname"]).strip()
        return res in water_resnames or res in ions

    ref_non_solvent = [dict(r) for r in ref_records if not _is_solvent_or_ion(r)]
    tgt_non_solvent = [r for r in tgt_records if not _is_solvent_or_ion(r)]
    tgt_solvent_ions = [dict(r) for r in tgt_records if _is_solvent_or_ion(r)]

    if len(ref_non_solvent) != len(tgt_non_solvent):
        print(
            "⚠ Could not restore non-solvent atom labeling after genion: "
            f"reference has {len(ref_non_solvent)} non-solvent atoms but target has {len(tgt_non_solvent)}."
        )
        return False

    merged = ref_non_solvent + tgt_solvent_ions
    for i, rec in enumerate(merged, start=1):
        rec["atomid"] = i
    _write_gro_records(str(target_gro), title, merged, box)
    return True


def _extract_molecules_block(top_path: Path) -> str:
    text = top_path.read_text()
    lines = text.splitlines()
    start = None
    for i, line in enumerate(lines):
        if line.strip().lower() == "[ molecules ]":
            start = i
            break
    if start is None:
        raise RuntimeError(f"[ molecules ] block not found in {top_path}")
    block = "\n".join(lines[start:]).rstrip() + "\n"
    return block


def _parse_molecules_entries(top_path: Path) -> list[tuple[str, int]]:
    text = top_path.read_text().splitlines()
    in_molecules = False
    entries: list[tuple[str, int]] = []
    for raw in text:
        s = raw.strip()
        if s.lower() == "[ molecules ]":
            in_molecules = True
            continue
        if not in_molecules:
            continue
        if s.startswith("["):
            break
        if not s or s.startswith(";"):
            continue
        parts = s.split()
        if len(parts) < 2:
            continue
        name = parts[0]
        try:
            count = int(float(parts[1]))
        except ValueError:
            continue
        entries.append((name, count))
    return entries


def _normalize_molecules_block(
    final_top: Path,
    base_top: Path | None = None,
    top_dir: Path | None = None,
) -> None:
    final_entries = _parse_molecules_entries(final_top)
    if not final_entries:
        return

    # Deduplicate by molname (keep latest count), while preserving first-seen order.
    counts: dict[str, int] = {}
    first_order: list[str] = []
    for name, count in final_entries:
        if name not in counts:
            first_order.append(name)
        counts[name] = int(count)

    base_order: list[str] = []
    base_counts: dict[str, int] = {}
    if base_top is not None and base_top.exists():
        for name, bcount in _parse_molecules_entries(base_top):
            if name not in base_order:
                base_order.append(name)
            # Keep first definition from base topology.
            if name not in base_counts:
                base_counts[name] = int(bcount)

    # Safety net: keep fundamental molecules from base topology even if
    # downstream tools accidentally dropped them from final_top.
    for name, bcount in base_counts.items():
        if name not in counts:
            counts[name] = int(bcount)
            first_order.append(name)

    ff_ions = _detect_ff_ion_moltypes(top_dir) if top_dir is not None else {"NA": "NA", "CL": "CL"}

    # Canonicalize sodium/chloride names to the active FF moltypes.
    ion_alias_pairs = [("NA", "NA+"), ("CL", "CL-")]
    for base_name, alt_name in ion_alias_pairs:
        preferred = ff_ions.get(base_name, base_name)
        alias = alt_name if preferred == base_name else base_name
        if alias in counts:
            counts[preferred] = int(counts.get(preferred, 0)) + int(counts[alias])
            del counts[alias]

    water_names = ["W", "SOL", "WF"]
    ion_names = [
        ff_ions.get("NA", "NA"),
        ff_ions.get("CL", "CL"),
        "K", "CA", "MG", "ZN", "LI", "RB", "CS", "BA", "SR", "F", "BR", "I",
    ]

    # Non-solvent molecules should keep base topology counts (DNA/protein/surface/linker/cofactor/substrate).
    for name, bcount in base_counts.items():
        if name not in water_names and name not in ion_names:
            counts[name] = int(bcount)

    ordered: list[str] = []

    # 1) Keep non-solvent/non-ion molecules in base topology order (e.g., biomolecule, surface, linker, substrate).
    for name in base_order:
        if name in counts and counts[name] > 0 and name not in water_names and name not in ion_names:
            ordered.append(name)

    # 2) Append any remaining non-solvent/non-ion molecules in first appearance order.
    for name in first_order:
        if name in counts and counts[name] > 0 and name not in water_names and name not in ion_names and name not in ordered:
            ordered.append(name)

    # 3) Waters.
    for name in water_names:
        if name in counts and counts[name] > 0:
            ordered.append(name)

    # 4) Ions.
    for name in ion_names:
        if name in counts and counts[name] > 0:
            ordered.append(name)

    # 5) Any unknown leftovers.
    for name in first_order:
        if name in counts and counts[name] > 0 and name not in ordered:
            ordered.append(name)

    block_lines = ["[ molecules ]"]
    for name in ordered:
        block_lines.append(f"{name} {counts[name]}")
    molecules_block = "\n".join(block_lines) + "\n"
    _replace_molecules_block(final_top, molecules_block)


def _replace_molecules_block(top_path: Path, molecules_block: str) -> None:
    text = top_path.read_text()
    lines = text.splitlines()
    start = None
    for i, line in enumerate(lines):
        if line.strip().lower() == "[ molecules ]":
            start = i
            break
    if start is None:
        top_path.write_text(text.rstrip() + "\n\n" + molecules_block)
        return
    head = "\n".join(lines[:start]).rstrip() + "\n\n"
    top_path.write_text(head + molecules_block)


def _sync_final_restrained_topology(top_dir: Path, final_top: Path) -> Path | None:
    base_res_top = top_dir / "system_res.top"
    if not base_res_top.exists():
        return None
    final_res_top = top_dir / "system_final_res.top"
    shutil.copy(base_res_top, final_res_top)
    molecules_block = _extract_molecules_block(final_top)
    _replace_molecules_block(final_res_top, molecules_block)
    return final_res_top


def _rebuild_merged_index(
    gmx_bin: str,
    gro_path: Path,
    top_dir: Path,
) -> Path:
    index_path = top_dir / "index.ndx"
    auto_index_path = top_dir / "_index_auto.ndx"

    custom_index_text = ""
    if index_path.exists():
        custom_index_text = index_path.read_text()

    _run_with_check(
        [
            gmx_bin, "make_ndx",
            "-f", str(gro_path),
            "-o", str(auto_index_path),
        ],
        cwd=top_dir,
        stdin_text="q\n",
    )

    auto_index_text = auto_index_path.read_text().rstrip()
    merged = auto_index_text + "\n"
    if custom_index_text.strip():
        merged += "\n; MartiniSurf custom index groups (appended)\n"
        merged += custom_index_text.strip() + "\n"
    index_path.write_text(merged)
    auto_index_path.unlink(missing_ok=True)
    return index_path


def _run_optional_solvation_ionization(args: argparse.Namespace, simdir: Path) -> None:
    if not args.solvate:
        return

    gmx_bin = _find_gmx_binary()
    if gmx_bin is None:
        raise RuntimeError(
            "Solvation/Ionization requested, but no GROMACS binary was found. "
            "Install GROMACS and make sure `gmx` is on PATH."
        )

    system_dir = simdir / "2_system"
    top_dir = simdir / "0_topology"

    input_gro = system_dir / "system.gro"
    if not input_gro.exists():
        input_gro = system_dir / "immobilized_system.gro"
    if not input_gro.exists():
        raise FileNotFoundError("Could not find system GRO for solvation.")

    if args.water_gro:
        water_gro = Path(args.water_gro).resolve()
    else:
        water_gro = system_dir / "water.gro"
        if not water_gro.exists():
            pkg_water = Path(__file__).resolve().parent / "system_templates" / "water.gro"
            if pkg_water.exists():
                water_gro = pkg_water
            else:
                raise FileNotFoundError("water.gro template not found for solvation.")

    base_top = top_dir / "system.top"
    if not base_top.exists():
        raise FileNotFoundError("system.top not found for solvation.")

    final_top = top_dir / "system_final.top"
    shutil.copy(base_top, final_top)

    solvated_gro = system_dir / "solvated_system.gro"
    _run_with_check([
        gmx_bin, "solvate",
        "-cp", str(input_gro),
        "-cs", str(water_gro),
        "-o", str(solvated_gro),
        "-p", str(final_top),
        "-radius", str(args.solvate_radius),
    ], cwd=top_dir)
    surface_itp = top_dir / "system_itp" / "surface.itp"
    surface_resname = _read_itp_moleculetype(surface_itp) or "SRF"
    water_resname = _read_gro_first_resname(str(water_gro)) or "W"
    water_resnames = {water_resname}
    if water_resname == "W":
        water_resnames.add("SOL")
    _remove_waters_below_surface(
        gro_path=solvated_gro,
        top_path=final_top,
        surface_resname=surface_resname,
        clearance_nm=args.solvate_surface_clearance,
        water_resnames=water_resnames,
        water_molname=water_resname,
    )

    final_gro = system_dir / "final_system.gro"
    final_alias_gro = system_dir / "system_final.gro"
    if not args.ionize:
        _normalize_molecules_block(final_top=final_top, base_top=base_top, top_dir=top_dir)
        shutil.copy(solvated_gro, final_gro)
        _normalize_uniform_atom_names_from_itp(top_dir=top_dir, top_path=final_top, gro_path=final_gro)
        shutil.copy(final_gro, final_alias_gro)
        merged_index = _rebuild_merged_index(gmx_bin=gmx_bin, gro_path=final_alias_gro, top_dir=top_dir)
        final_res_top = _sync_final_restrained_topology(top_dir, final_top)
        print(f"✔ Solvated system written: {final_gro}")
        print(f"✔ Alias final system written: {final_alias_gro}")
        print(f"✔ Final topology written: {final_top}")
        print(f"✔ Merged index written: {merged_index}")
        if final_res_top is not None:
            print(f"✔ Final restrained topology written: {final_res_top}")
        return

    ions_mdp = system_dir / "_ions_tmp.mdp"
    ions_tpr = system_dir / "_ions_tmp.tpr"
    _write_ions_mdp(ions_mdp)

    _run_with_check([
        gmx_bin, "grompp",
        "-f", str(ions_mdp),
        "-c", str(solvated_gro),
        "-p", str(final_top),
        "-o", str(ions_tpr),
        "-maxwarn", "2",
    ], cwd=top_dir)

    _run_genion_with_fallback(
        gmx_bin=gmx_bin,
        tpr_path=ions_tpr,
        out_gro=final_gro,
        top_path=final_top,
        salt_conc=args.salt_conc,
        cwd=top_dir,
    )
    _restore_non_solvent_from_reference(
        reference_gro=solvated_gro,
        target_gro=final_gro,
        water_resnames=water_resnames,
        ion_resnames={"NA", "CL"},
    )
    _normalize_uniform_atom_names_from_itp(top_dir=top_dir, top_path=final_top, gro_path=final_gro)
    _normalize_ion_atom_names_from_itp(top_dir=top_dir, gro_path=final_gro)
    _normalize_molecules_block(final_top=final_top, base_top=base_top, top_dir=top_dir)

    ions_mdp.unlink(missing_ok=True)
    ions_tpr.unlink(missing_ok=True)
    shutil.copy(final_gro, final_alias_gro)
    merged_index = _rebuild_merged_index(gmx_bin=gmx_bin, gro_path=final_alias_gro, top_dir=top_dir)
    final_res_top = _sync_final_restrained_topology(top_dir, final_top)
    print(f"✔ Ionized system written: {final_gro}")
    print(f"✔ Alias final system written: {final_alias_gro}")
    print(f"✔ Final topology written: {final_top}")
    print(f"✔ Merged index written: {merged_index}")
    if final_res_top is not None:
        print(f"✔ Final restrained topology written: {final_res_top}")

    # Keep topology includes centralized only in 0_topology/system_itp.
    accidental_itp_dir = system_dir / "system_itp"
    if accidental_itp_dir.exists():
        shutil.rmtree(accidental_itp_dir, ignore_errors=True)


def _run_optional_dna_water_freezing(args: argparse.Namespace, simdir: Path) -> None:
    if args.freeze_water_fraction <= 0:
        return

    from martinisurf.utils.freeze_water import apply_freeze_water_fraction

    top_dir = simdir / "0_topology"
    system_dir = simdir / "2_system"
    top_path = top_dir / "system_final.top"
    gro_path = system_dir / "final_system.gro"
    alias_path = system_dir / "system_final.gro"

    n_before, n_w, n_wf = apply_freeze_water_fraction(
        top_path=top_path,
        gro_path=gro_path,
        fraction=args.freeze_water_fraction,
        seed=args.freeze_water_seed,
        source_resname="W",
        target_resname="WF",
        alias_gro_path=alias_path,
    )
    _normalize_molecules_block(final_top=top_path, base_top=top_dir / "system.top", top_dir=top_dir)
    final_res_top = _sync_final_restrained_topology(top_dir, top_path)
    gmx_bin = _find_gmx_binary()
    if gmx_bin is not None:
        merged_index = _rebuild_merged_index(gmx_bin=gmx_bin, gro_path=alias_path, top_dir=top_dir)
        print(f"✔ Merged index written: {merged_index}")
    if final_res_top is not None:
        print(f"✔ Final restrained topology written: {final_res_top}")
    print(
        "✔ DNA water freezing applied "
        f"(W before={n_before}, W after={n_w}, WF={n_wf}, seed={args.freeze_water_seed})"
    )


def _run_final_topology_structure_validation(simdir: Path) -> None:
    top_dir = simdir / "0_topology"
    system_dir = simdir / "2_system"

    top_path = top_dir / "system_final_res.top"
    if not top_path.exists():
        top_path = top_dir / "system_final.top"
    gro_path = system_dir / "system_final.gro"
    if not gro_path.exists():
        gro_path = system_dir / "final_system.gro"

    if not top_path.exists() or not gro_path.exists():
        return

    _validate_named_molecule_atomnames(top_dir=top_dir, top_path=top_path, gro_path=gro_path)


def _parse_version_triplet(text: str) -> tuple[int, int, int] | None:
    match = re.search(r"(\d+)\.(\d+)\.(\d+)", text or "")
    if not match:
        return None
    return int(match.group(1)), int(match.group(2)), int(match.group(3))


def _version_leq(found: tuple[int, int, int], limit: tuple[int, int, int]) -> bool:
    return found <= limit


def _is_dssp_binary_compatible(binary_path: str) -> bool:
    try:
        res = subprocess.run(
            [binary_path, "--version"],
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError:
        return False

    version_text = f"{res.stdout}\n{res.stderr}"
    parsed = _parse_version_triplet(version_text)
    if parsed is None:
        # If version cannot be parsed, avoid hard-failing and let martinize2 retry logic handle it.
        print(f"⚠ Could not parse DSSP version for {binary_path}.")
        return True

    compatible = _version_leq(parsed, (3, 1, 4))
    if not compatible:
        print(
            f"⚠ DSSP at {binary_path} is version {parsed[0]}.{parsed[1]}.{parsed[2]}, "
            "but martinize2 is only compatible with DSSP <= 3.1.4."
        )
    return compatible


def _select_dssp_flags() -> list[str]:
    # Martinize2 recommendation: prefer mdtraj for secondary structure assignment.
    if importlib.util.find_spec("mdtraj") is not None:
        return ["-dssp"]

    # Fallback to binary DSSP only when compatible.
    dssp_env = os.environ.get("DSSP", "").strip()
    if dssp_env and _is_dssp_binary_compatible(dssp_env):
        return ["-dssp", dssp_env]

    system_mkdssp = shutil.which("mkdssp")
    if system_mkdssp and _is_dssp_binary_compatible(system_mkdssp):
        return ["-dssp", system_mkdssp]

    dssp_bin = Path(__file__).resolve().parent / "dssp" / "mkdssp"
    if dssp_bin.exists() and _is_dssp_binary_compatible(str(dssp_bin)):
        return ["-dssp", str(dssp_bin)]

    print("⚠ DSSP requested but no compatible setup found. Continuing without DSSP.")
    print("  Install `mdtraj` (recommended) or provide DSSP <= 3.1.4 via $DSSP.")
    return []


def _backup_existing_output_dir(simdir: Path) -> Path | None:
    if not simdir.exists():
        return None

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup = simdir.with_name(f"{simdir.name}_backup_{timestamp}")
    suffix = 1
    while backup.exists():
        backup = simdir.with_name(f"{simdir.name}_backup_{timestamp}_{suffix}")
        suffix += 1

    shutil.move(str(simdir), str(backup))
    print(f"ℹ Existing output moved to backup: {backup}")
    return backup


# ======================================================================
# MAIN
# ======================================================================

def main(argv=None):

    parser = build_parser()
    args = parser.parse_args(argv)
    _validate_args(parser, args)
    _print_config_summary(args)
    merge_groups = _normalize_merge_groups(args.merge)

    simdir = Path(args.outdir).resolve()
    _backup_existing_output_dir(simdir)

    (simdir / "0_topology" / "system_itp").mkdir(parents=True)
    (simdir / "1_mdp").mkdir()
    (simdir / "2_system").mkdir()

    active_itp_dir = simdir / "0_topology" / "system_itp"
    system_dir     = simdir / "2_system"

    complex_cfg: dict[str, Any] | None = None
    tmpdir: Path | None = None
    protein_go_model = bool(args.go)
    if args.complex_config:
        complex_cfg = _load_pre_cg_complex_config(Path(args.complex_config).resolve())
        mol = complex_cfg["protein_molname"]
        protein_go_model = bool(complex_cfg["include_go"] or args.go)

        complex_gro_path: Path = complex_cfg["complex_gro"]
        shutil.copy(complex_gro_path, system_dir / f"{mol}_cg.gro")

        protein_src_itp: Path = complex_cfg["protein_itp"]
        protein_dst_itp = active_itp_dir / f"{mol}.itp"
        shutil.copy(protein_src_itp, protein_dst_itp)

        cofactor_src_itp: Path = complex_cfg["cofactor_itp"]
        shutil.copy(cofactor_src_itp, active_itp_dir / cofactor_src_itp.name)

        if protein_go_model:
            if not complex_cfg["go_files"]:
                raise FileNotFoundError(
                    "Go model requested in pre-CG complex mode, but no go_* files were found. "
                    "Add them to input/ or adjust topology.go_files_glob."
                )
            for go_file in complex_cfg["go_files"]:
                shutil.copy(go_file, active_itp_dir / go_file.name)
    else:
        tmpdir = simdir / "_martinize_tmp"
        tmpdir.mkdir()

        mol = args.moltype if args.moltype else "DNA"

        system_cg_out = tmpdir / f"{mol}_cg.pdb"
        topfile_out   = tmpdir / f"{mol}_cg.top"

        # ===============================================================
        # 1) CLEAN INPUT
        # ===============================================================
        pdb_abs = load_clean_pdb(args.pdb, workdir=simdir)
        resolved_anchor_groups = _normalize_cli_residue_groups(args.anchor, pdb_abs, "--anchor")
        resolved_linker_groups = _normalize_cli_residue_groups(args.linker_group, pdb_abs, "--linker-group")

        # ===============================================================
        # 2) MARTINIZATION
        # ===============================================================

        if args.dna:
            print("🧬 DNA mode → using martinize-dna.py")

            dna_script = Path(__file__).resolve().parent / "utils" / "martinize-dna.py"
            from martinisurf.utils.use_python2 import find_python2
            python2 = find_python2()

            dna_input_gro = tmpdir / "dna_input.gro"
            pdb_to_gro(str(pdb_abs), str(dna_input_gro))

            martinize_cmd = [
                python2,
                str(dna_script),
                "-f", str(dna_input_gro),
                "-x", str(system_cg_out),
                "-o", str(topfile_out),
                "-dnatype", args.dnatype,
                "-p", args.p.capitalize(),
                "-pf", str(args.pf),
            ]
            for group in merge_groups:
                martinize_cmd += ["-merge", group]

            if args.elastic:
                martinize_cmd += ["-elastic", "-ef", str(args.ef)]

        else:
            print("🧬 Protein mode → using martinize2")

            martinize_cmd = [
                "martinize2",
                "-f", str(pdb_abs),
                "-x", str(system_cg_out),
                "-o", str(topfile_out),
                "-ff", args.ff,
                "-name", mol,
                "-maxwarn", str(args.maxwarn),
            ]
            for group in merge_groups:
                martinize_cmd += ["-merge", group]

            if args.p != "none":
                martinize_cmd += ["-p", args.p]

            martinize_cmd += ["-pf", str(args.pf)]

            if args.elastic:
                martinize_cmd += ["-elastic", "-ef", str(args.ef)]

            if args.go:
                martinize_cmd += ["-go"]
                if args.go_eps is not None:
                    martinize_cmd += ["-go-eps", str(args.go_eps)]
                if args.go_low is not None:
                    martinize_cmd += ["-go-low", str(args.go_low)]
                if args.go_up is not None:
                    martinize_cmd += ["-go-up", str(args.go_up)]

            if args.dssp:
                martinize_cmd += _select_dssp_flags()

        # martinize2 in some Colab/runtime setups fails when DSSP binary is present
        # but not functional. In that case retry once without DSSP.
        try:
            run(martinize_cmd, cwd=tmpdir)
        except RuntimeError:
            if (not args.dna) and args.dssp and "-dssp" in martinize_cmd:
                print("⚠ martinize2 failed with DSSP. Retrying without DSSP...")
                retry_cmd = martinize_cmd[:]
                dssp_idx = retry_cmd.index("-dssp")
                # Remove -dssp and optional binary path if present.
                del retry_cmd[dssp_idx]
                if dssp_idx < len(retry_cmd) and not retry_cmd[dssp_idx].startswith("-"):
                    del retry_cmd[dssp_idx]
                run(retry_cmd, cwd=tmpdir)
            else:
                raise

        # Move ITP files
        for f in tmpdir.glob("*.itp"):
            shutil.move(str(f), active_itp_dir / f.name)

        shutil.copy(system_cg_out, system_dir / f"{mol}_cg.pdb")
    if args.complex_config:
        resolved_anchor_groups = complex_cfg["anchor_groups"]
        resolved_linker_groups = []

    # ===============================================================
    # 3) SURFACE
    # ===============================================================

    surface_gro = system_dir / "surface.gro"

    if args.surface:

        # Copy GRO
        shutil.copy(args.surface, surface_gro)

        # 🔹 Copy surface.itp if exists
        possible_itp = Path(args.surface).with_suffix(".itp")
        if possible_itp.exists():
            shutil.copy(possible_itp, active_itp_dir / "surface.itp")
            print("✔ Copied surface.itp from provided surface file")
        else:
            fallback_itp = active_itp_dir / "surface.itp"
            fallback_resname = _read_gro_first_resname(str(surface_gro)) or "SRF"
            fallback_bead = _read_gro_first_atomname(str(surface_gro)) or args.surface_bead
            _write_minimal_surface_itp(
                itp_path=fallback_itp,
                resname=fallback_resname,
                bead=fallback_bead,
                charge=float(args.charge),
            )
            print(
                "⚠ No surface.itp found next to provided surface.gro. "
                f"Generated fallback topology: {fallback_itp}"
            )

    else:
        import martinisurf.surface_builder as sb

        sb.main([
            "--mode", args.surface_mode,
            "--lx", str(args.lx),
            "--ly", str(args.ly),
            "--dx", str(args.dx),
            "--bead", args.surface_bead,
            "--charge", str(args.charge),
            "--output", str(system_dir / "surface"),
        ])

        # 🔹 Move generated surface.itp into topology
        generated_itp = system_dir / "surface.itp"
        if generated_itp.exists():
            shutil.move(generated_itp, active_itp_dir / "surface.itp")
            print("✔ Moved generated surface.itp into topology")
        else:
            print("⚠ surface.itp not generated by surface_builder")

    # Keep topology includes centralized only in 0_topology/system_itp.
    accidental_itp_dir = system_dir / "system_itp"
    if accidental_itp_dir.exists():
        shutil.rmtree(accidental_itp_dir, ignore_errors=True)

    # ===============================================================
    # 4) ORIENTATION
    # ===============================================================

    import martinisurf.system_tethered as orient_mod

    system_gro = system_dir / f"{mol}_cg.gro"
    if not complex_cfg:
        pdb_to_gro(str(system_dir / f"{mol}_cg.pdb"), str(system_gro))

    orient_args = [
        "--surface", str(surface_gro),
        "--system", str(system_gro),
        "--out", str(system_dir / "immobilized_system.gro"),
    ]
    if args.dna:
        orient_args += ["--dna-mode"]
    balance_low_z = bool(args.balance_low_z or (complex_cfg and complex_cfg.get("balance_low_z")))
    balance_low_z_fraction = args.balance_low_z_fraction
    if balance_low_z_fraction is None and complex_cfg:
        balance_low_z_fraction = complex_cfg.get("balance_low_z_fraction")
    if balance_low_z:
        orient_args += ["--balance-low-z"]
    if balance_low_z_fraction is not None:
        orient_args += ["--balance-low-z-fraction", str(balance_low_z_fraction)]

    if args.linker:
        first_group_resid = None
        if resolved_linker_groups and len(resolved_linker_groups[0]) > 1:
            try:
                first_group_resid = int(resolved_linker_groups[0][1])
            except ValueError:
                first_group_resid = None

        linker_first = _read_gro_first_atomname(args.linker)
        linker_last = _read_gro_last_atomname(args.linker)
        # Protein side is the first linker bead by default; swap when inverted.
        linker_head = linker_last if args.invert_linker else linker_first
        linker_tail = linker_first if args.invert_linker else linker_last
        target_atom = _read_gro_atomname_for_resid(str(system_gro), first_group_resid) if first_group_resid else None
        surface_atom = _read_gro_first_atomname(str(surface_gro))

        prot_sigma_nm = _sigma_nm(
            is_dna=args.dna,
            ff_name=args.ff,
            class_a=_bead_size_class(linker_head),
            class_b=_bead_size_class(target_atom),
        )
        surf_sigma_nm = _sigma_nm(
            is_dna=args.dna,
            ff_name=args.ff,
            class_a=_bead_size_class(linker_tail),
            class_b=_bead_size_class(surface_atom or args.surface_bead),
        )

        # Auto distances in nm:
        # - DNA linker-to-biomolecule uses raw sigma (bonded coupling in topology).
        # - Protein linker-to-biomolecule keeps legacy sigma*1.2 pull distance.
        linker_prot_dist_nm = (
            args.linker_prot_dist
            if args.linker_prot_dist is not None
            else (prot_sigma_nm if args.dna else prot_sigma_nm * 1.2)
        )
        linker_surf_dist_nm = (
            args.linker_surf_dist
            if args.linker_surf_dist is not None
            else surf_sigma_nm * 1.2
        )
        linker_prot_dist_ang = linker_prot_dist_nm * 10.0
        linker_surf_dist_ang = linker_surf_dist_nm * 10.0

        print(
            f"ℹ Linker distances (nm): prot={linker_prot_dist_nm:.3f} "
            f"surf={linker_surf_dist_nm:.3f}"
        )

        orient_args += [
            "--linker-gro", str(args.linker),
            "--linker-prot-dist", str(linker_prot_dist_ang),
            "--linker-surf-dist", str(linker_surf_dist_ang),
        ]
        if args.invert_linker:
            orient_args += ["--invert-linker"]

        if args.surface_linkers > 0:
            orient_args += ["--surface-linkers", str(args.surface_linkers)]

        if resolved_linker_groups:
            for group in resolved_linker_groups:
                orient_args += ["--linker-group"] + [str(x) for x in group]

    elif args.anchor:
        orient_args += ["--dist", str(args.dist * 10.0)]
        for group in resolved_anchor_groups:
            orient_args += ["--anchor"] + [str(x) for x in group]
    elif complex_cfg:
        if not complex_cfg["anchor_groups"]:
            raise ValueError(
                "complex_config.yaml must define protein.anchor_groups or protein.orient_by_residues "
                "when not using --linker mode."
            )
        if complex_cfg.get("anchor_landmark_mode") != "residue":
            orient_args += ["--anchor-landmark-mode", str(complex_cfg["anchor_landmark_mode"])]
        if complex_cfg.get("cofactor_molname"):
            orient_args += ["--reference-exclude-resname", str(complex_cfg["cofactor_molname"])]
        orient_args += ["--min-reference-z-dist", str(args.dist * 10.0)]
        orient_args += ["--dist", str(args.dist * 10.0)]
        for group in complex_cfg["anchor_groups"]:
            orient_args += ["--anchor"] + [str(x) for x in group]

    orient_mod.main(orient_args)

    if args.substrate and args.substrate_count > 0:
        _append_random_substrates_to_gro(
            system_gro=str(system_dir / "immobilized_system.gro"),
            substrate_gro=str(args.substrate),
            substrate_count=args.substrate_count,
        )
        print(f"✔ Added {args.substrate_count} random substrate molecule(s) inside the box")

    # ===============================================================
    # 5) FINAL GMX SYSTEM (CRITICAL FIX)
    # ===============================================================

    import martinisurf.gromacs_inputs as gsm

    final_args = ["--outdir", str(simdir)]

    if not args.dna:
        final_args += ["--moltype", mol]
        if protein_go_model:
            final_args += ["--go-model"]

    if args.linker:
        final_args += ["--use-linker"]
        linker_resname = _read_gro_first_resname(args.linker)
        if linker_resname:
            final_args += ["--linker-resname", linker_resname]
        linker_size = _read_gro_atom_count(args.linker)
        if linker_size and linker_size > 0:
            final_args += ["--linker-size", str(linker_size)]
        final_args += ["--linker-pull-init-prot", str(linker_prot_dist_nm)]
        final_args += ["--linker-pull-init-surf", str(linker_surf_dist_nm)]

        linker_itp = Path(args.linker).with_suffix(".itp")
        if not linker_itp.exists():
            raise FileNotFoundError(
                f"Linker mode requires the matching .itp next to linker GRO: {linker_itp}"
            )
        shutil.copy(linker_itp, active_itp_dir / linker_itp.name)
        final_args += ["--linker-itp-name", linker_itp.name]
        print(f"✔ Copied {linker_itp.name} into topology")

    if args.substrate and args.substrate_count > 0:
        substrate_itp: Path | None = None
        if args.substrate_itp:
            substrate_itp = _resolve_sidecar_itp(args.substrate, args.substrate_itp, "Substrate")
        else:
            inferred_itp = Path(args.substrate).with_suffix(".itp")
            if inferred_itp.exists():
                substrate_itp = inferred_itp

        if substrate_itp is not None:
            shutil.copy(substrate_itp, active_itp_dir / substrate_itp.name)
            final_args += ["--substrate-itp-name", substrate_itp.name]
            print(f"✔ Copied {substrate_itp.name} into topology")
        else:
            substrate_moltype = _read_gro_first_resname(args.substrate)
            if not substrate_moltype:
                raise FileNotFoundError(
                    "Substrate ITP not found next to the GRO file, and the substrate resname "
                    "could not be inferred for Martini FF fallback."
                )
            final_args += ["--substrate-moltype", substrate_moltype]
            print(
                f"ℹ No substrate ITP found for {args.substrate}. "
                f"Will look for {substrate_moltype} in the Martini force-field includes."
            )
        final_args += ["--substrate-count", str(args.substrate_count)]

    if complex_cfg:
        cofactor_itp_path: Path = complex_cfg["cofactor_itp"]
        final_args += ["--cofactor-itp-name", cofactor_itp_path.name]
        final_args += ["--cofactor-count", str(complex_cfg["cofactor_count"])]

    # If classical anchor mode
    if (not args.ads_mode) and args.anchor:
        for group in resolved_anchor_groups:
            final_args += ["--anchor"] + [str(x) for x in group]

    # If linker mode → reuse linker groups as anchors
    elif (not args.ads_mode) and resolved_linker_groups:
        for group in resolved_linker_groups:
            final_args += ["--anchor"] + [str(x) for x in group]
    elif (not args.ads_mode) and complex_cfg:
        for group in complex_cfg["anchor_groups"]:
            final_args += ["--anchor"] + [str(x) for x in group]
    if args.ads_mode:
        final_args += ["--ads-mode"]

    gsm.main(final_args)
    _run_optional_solvation_ionization(args, simdir)
    _run_optional_dna_water_freezing(args, simdir)
    _run_final_topology_structure_validation(simdir)

    if tmpdir and tmpdir.exists():
        shutil.rmtree(tmpdir)

    print("\n=====================================")
    print("✔ MartiniSurf Finished Successfully")
    print("=====================================\n")


if __name__ == "__main__":
    main()
