# CLI Reference

MartiniSurf exposes one main command with module entrypoints.

## Compact package help

```bash
martinisurf -h
```

## Full package help

```bash
martinisurf --full-help
```

## Main pipeline help

```bash
martinisurf build -h
```

You can also call build flags directly without `build`:

```bash
martinisurf --pdb 1RJW --moltype Protein --surface-mode 4-1 --lx 15 --ly 15 --anchor B 8 10 11
```

For `--pdb` workflows, `--anchor` and `--linker-group` accept either legacy global residue ids or chain-based residue ids from the input PDB.
In `--complex-config` workflows, low-Z balancing is enabled by default (`protein.balance_low_z=true`, default fraction `0.2`).

## Adsorption mode (`--ads-mode`)

`--ads-mode` is an anchor-based adsorption workflow:
- Keeps classical anchor orientation behavior.
- Skips anchor pull/restraint topology generation.
- Generates MDP stages for `minimization`, `nvt`, `npt`, and `production` (no `deposition` stage).
- Is incompatible with linker mode (`--linker` / `--use-linker`).

## Module helps

```bash
martinisurf surface -h
martinisurf orient -h
martinisurf system -h
```

For complete flag-by-flag details, see {doc}`USER_GUIDE`.
