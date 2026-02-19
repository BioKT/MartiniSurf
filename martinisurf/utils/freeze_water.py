from __future__ import annotations

from pathlib import Path


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
    water_keys: list[tuple[int, str]] = []
    water_seen: set[tuple[int, str]] = set()
    for ln in atom_lines:
        if len(ln) < 44:
            continue
        resid = int(ln[0:5])
        resname = ln[5:10].strip()
        if resname == source_resname:
            key = (resid, resname)
            if key not in water_seen:
                water_seen.add(key)
                water_keys.append(key)

    n_source_total = len(water_keys)
    if n_source_total == 0:
        raise RuntimeError(f"No {source_resname} residues found in {gro_path}.")

    n_freeze = max(1, int(n_source_total * float(fraction)))
    # Deterministic freeze selection: convert the last W residues in GRO order.
    # This keeps unfrozen W first and frozen WF immediately after them.
    freeze_set = set(water_keys[-n_freeze:])

    new_atom_lines: list[str] = []
    for ln in atom_lines:
        if len(ln) < 44:
            new_atom_lines.append(ln)
            continue
        resid = int(ln[0:5])
        resname = ln[5:10].strip()
        if (resid, resname) in freeze_set:
            ln = ln[:5] + f"{target_resname:<5}" + f"{target_resname:>5}" + ln[15:]
        new_atom_lines.append(ln)

    new_lines = [lines[0], lines[1], *new_atom_lines, lines[-1]]
    gro_path.write_text("\n".join(new_lines) + "\n")
    if alias_gro_path is not None:
        alias_gro_path.write_text("\n".join(new_lines) + "\n")

    source_res = set()
    target_res = set()
    for ln in new_atom_lines:
        if len(ln) < 44:
            continue
        resid = int(ln[0:5])
        resname = ln[5:10].strip()
        if resname == source_resname:
            source_res.add((resid, resname))
        elif resname == target_resname:
            target_res.add((resid, resname))

    n_source_final = len(source_res)
    n_target_final = len(target_res)

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
