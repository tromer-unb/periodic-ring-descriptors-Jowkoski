"""Periodic graph construction utilities for the 3D descriptor."""

import networkx as nx
import numpy as np


def build_periodic_3d_graph(structure, cutoff=1.895, translations=(-1, 0, 1), trim_to_central=False):
    """Build a finite image graph that preserves 3D periodic connectivity.

    Nodes are tuples of site index and integer cell translations. The function
    expands the base-cell neighbor list into translated images and therefore
    avoids a quadratic search over all periodic copies.
    """
    lattice = structure.lattice.matrix
    nodes, positions, species, fracs = [], {}, {}, {}
    for tx in translations:
        for ty in translations:
            for tz in translations:
                shift = np.array([tx, ty, tz], dtype=float)
                for i, site in enumerate(structure):
                    node = (i, tx, ty, tz)
                    frac = np.asarray(site.frac_coords, dtype=float) + shift
                    nodes.append(node)
                    fracs[node] = frac
                    positions[node] = frac @ lattice
                    species[node] = str(site.specie)

    graph = nx.Graph()
    for node in nodes:
        graph.add_node(node, position=positions[node], species=species[node], frac=fracs[node])

    neighbors = structure.get_all_neighbors(cutoff)
    for a in nodes:
        i, tx, ty, tz = a
        for nb in neighbors[i]:
            j = int(nb.index)
            image = tuple(int(round(float(x))) for x in nb.image)
            b = (j, tx + image[0], ty + image[1], tz + image[2])
            if b in graph and a != b:
                graph.add_edge(a, b, distance=float(nb.nn_distance))

    if trim_to_central:
        central = {node for node in graph if node[1:] == (0, 0, 0)}
        keep = set(central)
        for node in central:
            keep.update(graph.neighbors(node))
        graph = graph.subgraph(keep).copy()
    return graph