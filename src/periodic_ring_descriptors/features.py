"""Utilities for converting ring-level output into flat structure features."""

from collections import Counter

import numpy as np


GEOMETRY_KEYS = (
    "area",
    "perimeter",
    "anisotropy",
    "mean_radius",
    "j_area",
    "j_perimeter",
    "j_anisotropy",
    "j_area_ratio",
    "j_perimeter_ratio",
)


def structure_feature_vector(result, ring_sizes=range(3, 21)):
    """Return a flat, model-ready dictionary from descriptor output.

    This convenience helper summarizes the common topology/geometry/Joukowsky
    blocks. It does not add dataset-specific defect chemistry or metadata used
    by the paper's benchmark models.
    """
    rings = result.get("rings", [])
    sizes = np.array([ring["ring_size"] for ring in rings], dtype=float)
    counts = Counter(int(x) for x in sizes)

    features = {
        "n_rings": int(len(rings)),
        "ring_size_mean": float(np.mean(sizes)) if len(sizes) else np.nan,
        "ring_size_std": float(np.std(sizes)) if len(sizes) else np.nan,
        "ring_size_min": float(np.min(sizes)) if len(sizes) else np.nan,
        "ring_size_max": float(np.max(sizes)) if len(sizes) else np.nan,
    }

    total = max(len(rings), 1)
    for size in ring_sizes:
        features[f"ring_{size}_fraction"] = counts.get(size, 0) / total

    for key in GEOMETRY_KEYS:
        values = np.array([ring[key] for ring in rings], dtype=float)
        features[f"{key}_mean"] = float(np.mean(values)) if len(values) else np.nan
        features[f"{key}_std"] = float(np.std(values)) if len(values) else np.nan

    for key in ("rings_per_atom", "ring_memberships_per_atom", "n_central_bonds"):
        if key in result:
            features[key] = result[key]
    return features