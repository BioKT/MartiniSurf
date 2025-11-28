def surface_color_from_atomname(name):
    if name.startswith("P"): return [1.0,0.55,0.0]
    if name.startswith("C"): return [0.5,0.5,0.5]
    if name.startswith("Q"):
        if name.endswith("n") or name.endswith("N"): return [1.0,0.0,0.0]
        if name.endswith("p") or name.endswith("P"): return [0.0,1.0,0.0]
    return [0.5,0.5,0.5]
