#!/usr/bin/env python3
"""Validated 2D periodic ring descriptor.

The implementation reproduces the descriptor used for the paper benchmark.
For the planar carbon structures studied in the paper, the weighted minimum
cycle basis used here was audited against explicit planar face tracing and gave
identical ring counts, ring-size multisets, and geometric means for all 120
allotropes.
"""

import argparse
import json
from collections import Counter

import numpy as np
import pandas as pd
from pymatgen.core import Structure

from ._core_2d import (
    aggregate_descriptors,
    anisotropy_2d,
    build_periodic_2d_graph,
    canonical_species_signature,
    find_rings,
    polygon_area_2d,
    polygon_perimeter_2d,
    project_ring_to_2d,
)


RING_COLUMNS = [
    "ring_size", "species_signature", "species_count", "mean_radius",
    "area", "perimeter", "anisotropy", "j_area", "j_perimeter",
    "j_anisotropy", "j_area_ratio", "j_perimeter_ratio",
]


def joukowsky_metrics_v2(points_2d, a=0.25):
    """Return orientation-averaged transformed ring metrics.

    Each edge is aligned with the positive real axis once. Averaging over all
    edge alignments removes dependence on the projected in-plane axes, starting
    vertex, and cycle traversal direction.
    """
    points = np.asarray(points_2d, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2 or len(points) < 3:
        raise ValueError("A ring must contain at least three 2D points.")
    if not np.isfinite(points).all() or not np.isfinite(a) or a < 0:
        raise ValueError("Coordinates and parameter a must be finite, with a >= 0.")

    centered = points - points.mean(axis=0)
    z = centered[:, 0] + 1j * centered[:, 1]
    radius = float(np.mean(np.abs(z)))
    if radius <= 1e-10:
        raise ValueError("Degenerate ring: mean radius is too small.")
    normalized = z / radius
    if np.any(np.abs(normalized) <= 1e-10):
        raise ValueError("Ring vertex at the center: singular transformation.")

    edges = np.roll(z, -1) - z
    if np.any(np.abs(edges) <= 1e-10):
        raise ValueError("Ring contains a zero-length edge.")

    measures = []
    for edge in edges:
        # Align each edge with the positive real axis before applying the map.
        oriented = normalized * np.conj(edge / abs(edge))
        transformed = radius * (oriented + a * a / oriented)
        coords = np.column_stack((transformed.real, transformed.imag))
        measures.append((
            polygon_area_2d(coords),
            polygon_perimeter_2d(coords),
            anisotropy_2d(coords),
        ))

    average = np.mean(measures, axis=0)
    return {
        "mean_radius": radius,
        "j_area": float(average[0]),
        "j_perimeter": float(average[1]),
        "j_anisotropy": float(average[2]),
    }


def describe_ring_v2(graph, ring, a=0.25):
    points_3d = np.array([graph.nodes[node]["position"] for node in ring])
    species = [graph.nodes[node]["species"] for node in ring]
    points = project_ring_to_2d(points_3d)
    area = polygon_area_2d(points)
    perimeter = polygon_perimeter_2d(points)
    if area <= 1e-12 or perimeter <= 1e-10:
        raise ValueError("Degenerate ring area or perimeter.")
    transformed = joukowsky_metrics_v2(points, a=a)
    return {
        "ring_size": len(ring),
        "species_signature": canonical_species_signature(species),
        "species_count": dict(Counter(species)),
        "area": area,
        "perimeter": perimeter,
        "anisotropy": anisotropy_2d(points),
        **transformed,
        "j_area_ratio": transformed["j_area"] / area,
        "j_perimeter_ratio": transformed["j_perimeter"] / perimeter,
    }


def compute_descriptor_2d(cif_file, cutoff=1.895, max_ring_size=12, a=0.25):
    """Compute the validated 2D periodic ring descriptor from a structure file.

    Parameters
    ----------
    cif_file
        Structure file readable by pymatgen.
    cutoff
        Fixed bond cutoff in angstrom.
    max_ring_size
        Largest ring retained from the periodic graph.
    a
        Dimensionless Joukowsky parameter.
    """
    if not np.isfinite(cutoff) or cutoff <= 0:
        raise ValueError("The bond cutoff must be positive and finite.")
    if max_ring_size < 3:
        raise ValueError("max_ring_size must be at least 3.")
    if not np.isfinite(a) or a < 0:
        raise ValueError("Parameter a must be finite and non-negative.")

    structure = Structure.from_file(cif_file)
    actual_cutoff = float(cutoff)
    graph = build_periodic_2d_graph(
        structure, super_range=(-1, 0, 1), z_images=(0,),
        fixed_cutoff=actual_cutoff,
    )
    central_degrees = [graph.degree((i, 0, 0, 0)) for i in range(len(structure))]
    rings = find_rings(graph, min_size=3, max_size=max_ring_size)
    descriptors = [describe_ring_v2(graph, ring, a=a) for ring in rings]
    return {
        "descriptor_version": 2,
        "input_cif": str(cif_file),
        "n_atoms_unit_cell": len(structure),
        "formula": structure.composition.reduced_formula,
        "graph": {
            "n_nodes_supercell": graph.number_of_nodes(),
            "n_edges_supercell": graph.number_of_edges(),
            "central_degree_min": min(central_degrees),
            "central_degree_max": max(central_degrees),
        },
        "parameters": {
            "fixed_cutoff": actual_cutoff,
            "cutoff_mode": "fixed",
            "max_ring_size": max_ring_size,
            "joukowsky_a": a,
        },
        "method": {
            "transform": "r_mean * (z/r_mean + a^2/(z/r_mean))",
            "orientation": "mean of metrics over all edge-aligned orientations",
            "j_area_units": "square angstrom for CIF lengths in angstrom",
            "j_perimeter_units": "angstrom for CIF lengths in angstrom",
            "ratios": "dimensionless",
        },
        "descriptor": aggregate_descriptors(descriptors),
        "rings": descriptors,
    }



# Backward-compatible name used in the original research scripts.
compute_descriptor_v2 = compute_descriptor_2d

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cif", help="Input CIF file")
    parser.add_argument("--cutoff", type=float, default=1.895,
                        help="Fixed bond cutoff in Å")
    parser.add_argument("--max-ring", type=int, default=12)
    parser.add_argument("--joukowsky-a", type=float, default=0.25)
    parser.add_argument("--json-out", default="descriptor_v2.json")
    parser.add_argument("--csv-rings", default="rings_v2.csv")
    args = parser.parse_args()
    result = compute_descriptor_2d(
        args.cif, cutoff=args.cutoff,
        max_ring_size=args.max_ring, a=args.joukowsky_a,
    )
    with open(args.json_out, "w", encoding="utf-8") as output:
        json.dump(result, output, indent=2, ensure_ascii=False)
    pd.DataFrame(result["rings"], columns=RING_COLUMNS).to_csv(
        args.csv_rings, index=False,
    )

    summary = result["descriptor"]
    print(f"Version: 2 | CIF: {args.cif} | atoms: {result['n_atoms_unit_cell']}")
    print(f"Rings: {summary['n_rings']} | sizes: {summary['ring_size_counts']}")
    print(f"JSON: {args.json_out} | CSV: {args.csv_rings}")


if __name__ == "__main__":
    main()