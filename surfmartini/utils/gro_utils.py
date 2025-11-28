import numpy as np

def load_gro(filename):
    coords, atoms = [], []
    lines = open(filename).readlines()[2:-1]
    for line in lines:
        try:
            resid = int(line[0:5])
            resname = line[5:10].strip()
            atomname = line[10:15].strip()
            atomid = int(line[15:20])
            x = float(line[20:28]) * 10
            y = float(line[28:36]) * 10
            z = float(line[36:44]) * 10
            atoms.append((resid,resname,atomname,atomid))
            coords.append([x,y,z])
        except:
            continue
    return np.array(coords), atoms
