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

## Module helps

```bash
martinisurf surface -h
martinisurf orient -h
martinisurf system -h
```

For complete flag-by-flag details, see {doc}`USER_GUIDE`.
