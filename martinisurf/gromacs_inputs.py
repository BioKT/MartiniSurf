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


# ======================================================================
# UTILITIES
# ======================================================================

def ensure_dir(path: str | Path) -> None:
    Path(path).mkdir(parents=True, exist_ok=True)


def write_list(values: List[int], fh, chunk: int = 15) -> None:
    for i in range(0, len(values), chunk):
        fh.write(" ".join(str(v) for v in values[i: i + chunk]) + "\n")


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
    linker_itp_name: str = "linker.itp",
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
) -> None:
    """Create simple master topology files expected by legacy workflows/tests."""
    if is_dna:
        # DNA FF files are mutually exclusive (both define [ defaults ]).
        # Include exactly one base FF, preferring v2.1 non-polarizable.
        dna_ff = None
        for candidate in ("martini_v2.1-dna.itp", "martini_v2.1P-dna.itp"):
            if (dst_itp_dir / candidate).exists():
                dna_ff = candidate
                break
        forcefield_itps = [dna_ff, "martini_v2.0_ions.itp"] if dna_ff else ["martini_v2.0_ions.itp"]
    else:
        forcefield_itps = [
            "martini_v3.0.0.itp",
            "martini_v3.0.0_solvents_v1.itp",
            "martini_v3.0.0_ions_v1.itp",
        ]

    present_forcefields = [itp for itp in forcefield_itps if itp and (dst_itp_dir / itp).exists()]

    def _include_block(molecule_itp_name: str) -> str:
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
            lines.append(f'#include "system_itp/{linker_itp_name}"')
        if cofactor_itp_name and cofactor_count > 0:
            lines.append(f'#include "system_itp/{cofactor_itp_name}"')
        if substrate_itp_name and substrate_count > 0:
            lines.append(f'#include "system_itp/{substrate_itp_name}"')
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
            _include_block(mol_itp_name)
            + "\n\n[ system ]\nMartiniSurf system\n\n[ molecules ]\n"
            + f"{moltype} 1\n{cofactor_line}{linker_line}{surface_moltype} {surface_count}\n{substrate_line}"
        )

    with open(system_res_top, "w") as fh:
        fh.write(
            _include_block(anchor_itp_name)
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


def _infer_surface_pull_group(text: str, is_dna: bool) -> str:
    m = re.search(r"^\s*pull-group2-name\s*=\s*([A-Za-z0-9_+-]+)", text, flags=re.MULTILINE)
    if m:
        return m.group(1)
    return "SRF" if is_dna else "GRA"


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


def _build_pull_block(anchor_count: int, surface_group: str, k: float = 1000.0, init: float = 0.8) -> str:
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
        lines.append(f"pull-coord{i}-init         = {init:.1f}")
        lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def _build_linker_pull_block(
    linker_count: int,
    surface_group: str,
    k: float = 1000.0,
    init_prot: float = 0.8,
    init_surf: float = 0.8,
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
    lines.append(f"pull-ncoords             = {2 * linker_count}")

    coord_id = 1
    for i in range(linker_count):
        selected_gid = 3 * i + 1
        linker_head_gid = 3 * i + 2
        linker_tail_gid = 3 * i + 3

        # coord: biomolecule anchor <-> linker head (near biomolecule)
        lines.append(f"pull-coord{coord_id}-geometry     = distance")
        lines.append(f"pull-coord{coord_id}-groups       = {selected_gid} {linker_head_gid}")
        lines.append(f"pull-coord{coord_id}-type         = umbrella")
        lines.append(f"pull-coord{coord_id}-k            = {k:.1f}")
        lines.append(f"pull-coord{coord_id}-rate         = 0")
        lines.append(f"pull-coord{coord_id}-dim          = N N Y")
        lines.append(f"pull-coord{coord_id}-init         = {init_prot:.3f}")
        lines.append("")
        coord_id += 1

        # coord: linker tail (near surface) <-> surface
        lines.append(f"pull-coord{coord_id}-geometry     = distance")
        lines.append(f"pull-coord{coord_id}-groups       = {linker_tail_gid} {surface_gid}")
        lines.append(f"pull-coord{coord_id}-type         = umbrella")
        lines.append(f"pull-coord{coord_id}-k            = {k:.1f}")
        lines.append(f"pull-coord{coord_id}-rate         = 0")
        lines.append(f"pull-coord{coord_id}-dim          = N N Y")
        lines.append(f"pull-coord{coord_id}-init         = {init_surf:.3f}")
        lines.append("")
        coord_id += 1

    return "\n".join(lines).rstrip() + "\n"


def write_custom_mdp(
    src: Path,
    dst: Path,
    anchor_count: int,
    is_dna: bool,
    linker_pull: bool = False,
    linker_count: int = 1,
    linker_pull_init_prot: float = 0.8,
    linker_pull_init_surf: float = 0.8,
) -> None:
    text = src.read_text()
    has_pull = bool(re.search(r"^\s*pull\s*=", text, flags=re.MULTILINE))
    surface_group = _infer_surface_pull_group(text, is_dna=is_dna)
    clean = _strip_pull_section(text)

    if has_pull:
        if linker_pull:
            clean = clean.rstrip() + "\n" + _build_linker_pull_block(
                linker_count=linker_count,
                surface_group=surface_group,
                init_prot=linker_pull_init_prot,
                init_surf=linker_pull_init_surf,
            )
        elif anchor_count > 0:
            clean = clean.rstrip() + "\n" + _build_pull_block(anchor_count, surface_group)

    dst.write_text(clean)


def _validate_cli_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    if args.use_linker and args.linker_size is not None and args.linker_size <= 0:
        parser.error("--linker-size must be > 0.")
    if args.cofactor_count < 0:
        parser.error("--cofactor-count must be >= 0.")
    if args.cofactor_count > 0 and not args.cofactor_itp_name:
        parser.error("--cofactor-count requires --cofactor-itp-name.")
    if args.substrate_count < 0:
        parser.error("--substrate-count must be >= 0.")
    if args.substrate_count > 0 and not args.substrate_itp_name:
        parser.error("--substrate-count requires --substrate-itp-name.")


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
    parser.add_argument("--linker-pull-init-prot", type=float, default=0.8, help="Pull init (nm) for biomolecule↔linker.")
    parser.add_argument("--linker-pull-init-surf", type=float, default=0.8, help="Pull init (nm) for linker↔surface.")
    parser.add_argument("--go-model", action="store_true", help="Add GO_VIRT define to generated protein topologies.")
    parser.add_argument("--cofactor-itp-name", help="Cofactor topology filename in system_itp.")
    parser.add_argument("--cofactor-count", type=int, default=0, help="Number of cofactor molecules in the system.")
    parser.add_argument("--substrate-itp-name", help="Substrate topology filename in system_itp.")
    parser.add_argument("--substrate-count", type=int, default=0, help="Number of substrate molecules in the system.")

    args = parser.parse_args(argv) if argv else parser.parse_args()
    _validate_cli_args(parser, args)

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

    # ===============================================================
    # BUILD ANCHOR GROUPS
    # ===============================================================

    anchor_atoms: Dict[int, List[int]] = {}
    pull_anchor_atoms: Dict[int, List[int]] = {}
    linker_pull_enabled = False
    linker_atom_ids: set[int] = set()

    if args.use_linker and args.linker_resname:
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
            linker_pull_enabled = True
            print(f"  → Anchor_2 (linker head): atom {linker_head}")
            print(f"  → Anchor_3 (linker tail): atom {linker_tail}")

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
    water_template = pkg_dir / "system_templates" / "water.gro"
    if water_template.exists():
        shutil.copy(water_template, sys_dir / "water.gro")
        print("✔ Copied water.gro template into 2_system")
    else:
        print(f"ℹ water.gro template not found at {water_template}. Skipping water.gro copy.")

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
        required = ["martini_v2.0_ions.itp"]
        # Copy only one mutually-exclusive DNA FF file.
        for candidate in ("martini_v2.1-dna.itp", "martini_v2.1P-dna.itp"):
            if (src_itp_dir / candidate).exists():
                required.insert(0, candidate)
                break
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

    if args.use_linker and not (dst_itp_dir / args.linker_itp_name).exists():
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

    linker_total = 0
    linker_moltype_name: str | None = None
    if args.use_linker and args.linker_resname:
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
    substrate_moltype_name: str | None = None
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
    # Active_anchor.itp (UNCHANGED LOGIC)
    # ===============================================================
    active_anchor = dst_itp_dir / f"{mol_itp.stem}_anchor.itp"

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

    for gid, atoms in anchor_atoms.items():
        for atom in atoms:
            new_lines.append(f"{atom} 1 1000 1000 0\n")

    new_lines.append("#endif\n")

    with open(active_anchor, "w") as fout:
        fout.write("".join(new_lines))

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

        # copy MDPs
        for src_name, dst_name in mdp_files.items():
            src = mdp_pkg / src_name
            if src.exists():
                write_custom_mdp(
                    src=src,
                    dst=mdp_dir / dst_name,
                    anchor_count=len(pull_anchor_atoms),
                    is_dna=is_dna,
                    linker_pull=linker_pull_enabled,
                    linker_count=(len(pull_anchor_atoms) // 3) if linker_pull_enabled else 1,
                    linker_pull_init_prot=args.linker_pull_init_prot,
                    linker_pull_init_surf=args.linker_pull_init_surf,
                )
                print(f"  ✔ {src_name} → {dst_name}")
            else:
                print(f"⚠ Missing MDP template: {src_name}")

    anchor_moltype = mol_itp.stem if is_dna else moltype
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

    write_top_files(
        topo_dir=topo_dir,
        dst_itp_dir=dst_itp_dir,
        moltype=moltype,
        mol_itp_name=mol_itp.name,
        anchor_itp_name=f"{anchor_moltype}_anchor.itp",
        is_dna=is_dna,
        use_linker=args.use_linker,
        go_model=args.go_model,
        linker_itp_name=args.linker_itp_name,
        linker_moltype=linker_moltype_name,
        linker_count=linker_total,
        cofactor_itp_name=args.cofactor_itp_name if include_cofactor_itp else None,
        cofactor_moltype=cofactor_moltype_name,
        cofactor_count=args.cofactor_count,
        substrate_itp_name=args.substrate_itp_name if include_substrate_itp else None,
        substrate_moltype=substrate_moltype_name,
        substrate_count=args.substrate_count,
        surface_moltype=surface_moltype,
        surface_count=surface_count,
    )

    # ===============================================================
    # INDEX (UNCHANGED)
    # ===============================================================
    if anchor_atoms:
        with open(topo_dir / "index.ndx", "w") as ndx:
            groups_for_index = pull_anchor_atoms if pull_anchor_atoms else anchor_atoms
            for gid, atoms in groups_for_index.items():
                ndx.write(f"\n[ Anchor_{gid} ]\n")
                write_list(atoms, ndx)

    if args.go_model:
        print("✔ GoMartini system built")
    else:
        print("✔ MartiniSurf system built")


if __name__ == "__main__":
    main()
