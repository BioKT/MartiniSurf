from __future__ import annotations

import random
from pathlib import Path


def _parse_box_lengths(box_line: str) -> tuple[float, float, float] | None:
    fields = box_line.split()
    if len(fields) < 3:
        return None
    try:
        lengths = tuple(float(v) for v in fields[:3])
    except ValueError:
        return None
    if any(v <= 0 for v in lengths):
        return None
    return lengths


def _distance_sq(
    xyz_a: tuple[float, float, float],
    xyz_b: tuple[float, float, float],
    box_lengths: tuple[float, float, float] | None,
) -> float:
    dx = xyz_a[0] - xyz_b[0]
    dy = xyz_a[1] - xyz_b[1]
    dz = xyz_a[2] - xyz_b[2]

    if box_lengths is not None:
        lx, ly, lz = box_lengths
        dx -= lx * round(dx / lx)
        dy -= ly * round(dy / ly)
        dz -= lz * round(dz / lz)

    return dx * dx + dy * dy + dz * dz


def _rename_water_block(block_lines: list[str], target_resname: str) -> list[str]:
    renamed: list[str] = []
    for ln in block_lines:
        if len(ln) < 20:
            renamed.append(ln)
            continue
        renamed.append(ln[:5] + f"{target_resname:<5}" + f"{target_resname:>5}" + ln[15:])
    return renamed


def apply_freeze_water_fraction(
    top_path: Path,
    gro_path: Path,
    fraction: float,
    seed: int = 42,
    source_resname: str = "W",
    target_resname: str = "WF",
    alias_gro_path: Path | None = None,
) -> tuple[int, int, int]:
    def _fmt_molecule_line(name: str, count: int) -> str:
        # Match the fixed-width style used in MartiniSurf [ molecules ] blocks.
        return f"{name:<16}{int(count)}"

    if fraction <= 0:
        return 0, 0, 0

    if not top_path.exists():
        raise FileNotFoundError(f"Missing topology: {top_path}")
    if not gro_path.exists():
        raise FileNotFoundError(f"Missing GRO: {gro_path}")

    lines = gro_path.read_text().splitlines()
    if len(lines) < 3:
        raise RuntimeError(f"Invalid GRO format: {gro_path}")

    atom_lines = lines[2:-1]
    box_lengths = _parse_box_lengths(lines[-1])

    blocks: list[dict[str, object]] = []
    for ln in atom_lines:
        if len(ln) < 44:
            continue
        resid = int(ln[0:5])
        resname = ln[5:10].strip()
        xyz = (
            float(ln[20:28]),
            float(ln[28:36]),
            float(ln[36:44]),
        )

        if blocks and blocks[-1]["resid"] == resid and blocks[-1]["resname"] == resname:
            blocks[-1]["lines"].append(ln)
            blocks[-1]["coords"].append(xyz)
            continue

        blocks.append(
            {
                "id": len(blocks),
                "resid": resid,
                "resname": resname,
                "lines": [ln],
                "coords": [xyz],
            }
        )

    source_blocks = [b for b in blocks if b["resname"] == source_resname]
    target_blocks = [b for b in blocks if b["resname"] == target_resname]
    n_source_total = len(source_blocks)
    if n_source_total == 0:
        raise RuntimeError(f"No {source_resname} residues found in {gro_path}.")

    n_freeze = min(n_source_total, max(1, int(n_source_total * float(fraction))))

    source_centroids: list[tuple[float, float, float]] = []
    for block in source_blocks:
        coords = block["coords"]
        n_coords = len(coords)
        source_centroids.append(
            (
                sum(x for x, _, _ in coords) / n_coords,
                sum(y for _, y, _ in coords) / n_coords,
                sum(z for _, _, z in coords) / n_coords,
            )
        )

    rng = random.Random(seed)
    priorities = [rng.random() for _ in source_blocks]
    first_idx = rng.randrange(len(source_blocks))
    frozen_idx = [first_idx]
    frozen_idx_set = {first_idx}
    min_dist_sq = [float("inf")] * len(source_blocks)

    while len(frozen_idx) < n_freeze:
        last_idx = frozen_idx[-1]
        last_xyz = source_centroids[last_idx]
        for idx, xyz in enumerate(source_centroids):
            if idx in frozen_idx_set:
                continue
            dist_sq = _distance_sq(xyz, last_xyz, box_lengths)
            if dist_sq < min_dist_sq[idx]:
                min_dist_sq[idx] = dist_sq

        next_idx = max(
            (idx for idx in range(len(source_blocks)) if idx not in frozen_idx_set),
            key=lambda idx: (min_dist_sq[idx], priorities[idx], -idx),
        )
        frozen_idx.append(next_idx)
        frozen_idx_set.add(next_idx)

    unfrozen_source_blocks: list[list[str]] = []
    frozen_source_blocks: list[list[str]] = []
    for idx, block in enumerate(source_blocks):
        lines_for_block = block["lines"]
        if idx in frozen_idx_set:
            frozen_source_blocks.append(_rename_water_block(lines_for_block, target_resname))
        else:
            unfrozen_source_blocks.append(list(lines_for_block))

    final_target_blocks = [list(block["lines"]) for block in target_blocks] + frozen_source_blocks

    ordered_blocks: list[list[str]] = []
    inserted_water_block = False
    water_resnames = {source_resname, target_resname}
    for block in blocks:
        if block["resname"] in water_resnames:
            if not inserted_water_block:
                ordered_blocks.extend(unfrozen_source_blocks)
                ordered_blocks.extend(final_target_blocks)
                inserted_water_block = True
            continue
        ordered_blocks.append(list(block["lines"]))

    if not inserted_water_block:
        ordered_blocks.extend(unfrozen_source_blocks)
        ordered_blocks.extend(final_target_blocks)

    new_atom_lines: list[str] = []
    atom_id = 1
    for block_lines in ordered_blocks:
        for ln in block_lines:
            if len(ln) >= 20:
                ln = ln[:15] + f"{atom_id:5d}" + ln[20:]
            new_atom_lines.append(ln)
            atom_id += 1

    new_lines = [lines[0], lines[1], *new_atom_lines, lines[-1]]
    gro_path.write_text("\n".join(new_lines) + "\n")
    if alias_gro_path is not None:
        alias_gro_path.write_text("\n".join(new_lines) + "\n")

    n_source_final = len(unfrozen_source_blocks)
    n_target_final = len(final_target_blocks)

    text = top_path.read_text().splitlines()
    out: list[str] = []
    in_molecules = False
    seen_source = False
    seen_target = False
    inserted_target_after_source = False
    for raw in text:
        s = raw.strip()
        if s.lower() == "[ molecules ]":
            in_molecules = True
            out.append(raw)
            continue

        if in_molecules and s.startswith("["):
            if not seen_source:
                out.append(_fmt_molecule_line(source_resname, n_source_final))
                seen_source = True
            if n_target_final > 0 and not seen_target:
                out.append(_fmt_molecule_line(target_resname, n_target_final))
                seen_target = True
            in_molecules = False
            out.append(raw)
            continue

        if in_molecules and s and not s.startswith(";"):
            name = s.split()[0]
            if name == source_resname:
                out.append(_fmt_molecule_line(source_resname, n_source_final))
                seen_source = True
                if n_target_final > 0:
                    out.append(_fmt_molecule_line(target_resname, n_target_final))
                    seen_target = True
                    inserted_target_after_source = True
                continue
            if name == target_resname:
                if inserted_target_after_source:
                    continue
                out.append(_fmt_molecule_line(target_resname, n_target_final))
                seen_target = True
                continue

        out.append(raw)

    if in_molecules:
        if not seen_source:
            out.append(_fmt_molecule_line(source_resname, n_source_final))
            seen_source = True
        if n_target_final > 0 and not seen_target:
            out.append(_fmt_molecule_line(target_resname, n_target_final))

    top_path.write_text("\n".join(out).rstrip() + "\n")

    return n_source_total, n_source_final, n_target_final
