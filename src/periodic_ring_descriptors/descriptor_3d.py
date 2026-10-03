#!/usr/bin/env python3
"""Periodic 3D Joukowsky descriptor using cell-invariant edge-shortest rings.

The original experimental 3D implementation uses a finite cycle basis whose
selected cycles depend on the crystallographic cell. This version defines a
ring through each bond as the bond plus every shortest alternate path between
its endpoints in a periodic image graph. Translational copies are
canonicalized before aggregation.
"""
from collections import Counter
from pathlib import Path
import argparse
import json
import networkx as nx
import numpy as np
from pymatgen.core import Structure

from ._graph_3d import build_periodic_3d_graph
from .descriptor_2d import describe_ring_v2

DESCRIPTOR_VERSION = "3d-edge-shortest-v1"

def canonical_cycle_key(cycle):
    """Return a translation-, start-, and direction-invariant cycle key."""
    representations = []
    for sequence in (cycle, list(reversed(cycle))):
        for start in range(len(sequence)):
            rotated = sequence[start:] + sequence[:start]
            shift = np.array(rotated[0][1:], dtype=int)
            representation = tuple(
                (node[0], *(np.array(node[1:], dtype=int) - shift).tolist())
                for node in rotated
            )
            representations.append(representation)
    return min(representations)

def central_periodic_edges(graph, tol=1e-8):
    """Choose one periodic representative of each bond by midpoint position."""
    edges = []
    for u, v in graph.edges():
        fu = np.asarray(graph.nodes[u]["frac"], dtype=float)
        fv = np.asarray(graph.nodes[v]["frac"], dtype=float)
        midpoint = 0.5 * (fu + fv)
        if np.all(midpoint >= -tol) and np.all(midpoint < 1.0 - tol):
            edges.append((u, v))
    return edges

def find_edge_shortest_rings(graph, min_size=3, max_size=20):
    """Enumerate bond-shortest rings and remove translational duplicates."""
    rings = {}
    representative_edges = central_periodic_edges(graph)
    for u, v in representative_edges:
        attributes = dict(graph.get_edge_data(u, v))
        graph.remove_edge(u, v)
        try:
            alternate_length = nx.shortest_path_length(graph, u, v)
            ring_size = alternate_length + 1
            if min_size <= ring_size <= max_size:
                for path in nx.all_shortest_paths(graph, u, v):
                    cycle = list(path)
                    rings.setdefault(canonical_cycle_key(cycle), cycle)
        except nx.NetworkXNoPath:
            pass
        finally:
            graph.add_edge(u, v, **attributes)
    return list(rings.values()), len(representative_edges)

def compute_descriptor_3d_structure(structure, cutoff=1.895, max_ring_size=20,
                        joukowsky_a=0.25, image_range=1):
    """Compute the 3D descriptor for an already loaded pymatgen Structure."""
    translations = tuple(range(-image_range, image_range + 1))
    graph = build_periodic_3d_graph(
        structure, cutoff=cutoff, translations=translations,
        trim_to_central=False,
    )
    rings, n_bonds = find_edge_shortest_rings(
        graph, min_size=3, max_size=max_ring_size
    )
    descriptors = [
        describe_ring_v2(graph, ring, a=joukowsky_a) for ring in rings
    ]
    degrees = [graph.degree((i, 0, 0, 0)) for i in range(len(structure))]
    size_counts = Counter(item["ring_size"] for item in descriptors)
    memberships = sum(item["ring_size"] for item in descriptors)
    return {
        "descriptor_version": DESCRIPTOR_VERSION,
        "n_atoms": len(structure),
        "n_central_bonds": n_bonds,
        "degree_min": min(degrees),
        "degree_max": max(degrees),
        "central_degrees": degrees,
        "n_rings": len(descriptors),
        "ring_size_counts": dict(sorted(size_counts.items())),
        "rings_per_atom": len(descriptors) / len(structure),
        "ring_memberships_per_atom": memberships / len(structure),
        "ring_method": "bond plus every shortest alternate path",
        "parameters": {
            "cutoff_angstrom": float(cutoff),
            "max_ring_size": int(max_ring_size),
            "joukowsky_a": float(joukowsky_a),
            "image_range": int(image_range),
        },
        "rings": descriptors,
    }
def compute_descriptor_3d(cif, cutoff=1.895, max_ring_size=20,
              joukowsky_a=0.25, image_range=1):
    structure = Structure.from_file(cif)
    result = compute_descriptor_3d_structure(
        structure, cutoff=cutoff, max_ring_size=max_ring_size,
        joukowsky_a=joukowsky_a, image_range=image_range,
    )
    result["input_cif"] = str(cif)
    return result


# Backward-compatible names used in the original research scripts.
calculate_structure = compute_descriptor_3d_structure
calculate = compute_descriptor_3d

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cif")
    parser.add_argument("--cutoff", type=float, default=1.895)
    parser.add_argument("--max-ring", type=int, default=20)
    parser.add_argument("--joukowsky-a", type=float, default=0.25)
    parser.add_argument("--image-range", type=int, default=1)
    parser.add_argument("--json-out", default="descriptor_3d_v2.json")
    args = parser.parse_args()
    result = compute_descriptor_3d(
        args.cif, cutoff=args.cutoff, max_ring_size=args.max_ring,
        joukowsky_a=args.joukowsky_a, image_range=args.image_range,
    )
    Path(args.json_out).write_text(json.dumps(result, indent=2))
    print(
        f"{DESCRIPTOR_VERSION}: atoms={result['n_atoms']} "
        f"rings={result['n_rings']} sizes={result['ring_size_counts']}"
    )

if __name__ == "__main__":
    main()