import numpy as np


def compute_centroid(coords):
    return np.mean(coords, axis=0)


def pca_normal(coords):
    m = coords - coords.mean(axis=0)
    _, _, vh = np.linalg.svd(m)
    return vh[-1]
