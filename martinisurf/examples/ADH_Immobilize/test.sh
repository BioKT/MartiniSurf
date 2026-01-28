martinisurf --pdb 1RJW.pdb --moltype Protein --lx 12 --ly 12 --surface-bead P4 --anchor 1 8 11 --anchor 2 1025 1028 --dist 8 --merge A,B,C,D --maxwarn 1 --go --moltype Active
### Using surface previous generate

martinisurf --pdb 1RJW.pdb --go --moltype Protein --surface surface.gro --anchor 1 8 10 11 --anchor 2 1025 1027 1028 --dist 10 --merge A,B,C,D
