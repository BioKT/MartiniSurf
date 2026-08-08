from __future__ import annotations

from pathlib import Path
import os
from shutil import which

from .command_builder import BuildConfig


def validate_config(config: BuildConfig) -> list[str]:
    errors: list[str] = []

    if not str(config.input_path).strip():
        errors.append("Provide a protein structure file, RCSB ID, or UniProt ID.")

    if not config.surface_path and (config.lx <= 0 or config.ly <= 0):
        errors.append("Generated surfaces need positive X and Y dimensions.")

    if config.orientation_mode == "Linker":
        if not config.linker_path:
            errors.append("Linker mode needs a linker .gro file.")
        if not config.linker_groups:
            errors.append("Linker mode needs at least one linker group.")
    else:
        if not config.anchors:
            errors.append("Anchor or adsorption mode needs at least one anchor group.")

    if config.ionize and not config.solvate:
        errors.append("Ionization requires solvation.")

    if config.substrate_count > 0 and not config.substrate_path:
        errors.append("Substrate count is set, but no substrate .gro file was provided.")

    for path in [config.surface_path, config.linker_path, config.substrate_path, config.substrate_itp_path]:
        if isinstance(path, Path) and not path.exists():
            errors.append(f"Missing file: {path}")

    return errors


def _which_with_extra_paths(executable: str, extra_tool_dirs: list[Path] | None = None) -> str | None:
    path_entries = [str(path) for path in extra_tool_dirs or [] if path.exists()]
    path_entries.append(os.environ.get("PATH", ""))
    return which(executable, path=os.pathsep.join(path_entries))


def check_external_tools(
    config: BuildConfig,
    extra_tool_dirs: list[Path] | None = None,
) -> list[str]:
    warnings: list[str] = []

    if _which_with_extra_paths("martinize2", extra_tool_dirs) is None:
        warnings.append(
            "martinize2 is not available in PATH. Protein coarse-graining will fail until it is installed "
            "or the app is launched from an environment where martinize2 is active."
        )

    if (config.solvate or config.ionize) and _which_with_extra_paths("gmx", extra_tool_dirs) is None:
        warnings.append(
            "GROMACS (`gmx`) is not available in PATH. Solvation or ionization needs GROMACS."
        )

    return warnings
