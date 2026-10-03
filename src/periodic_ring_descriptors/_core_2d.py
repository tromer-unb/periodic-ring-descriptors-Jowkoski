#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Joukowsky Ring Descriptor para estruturas 2D tipo grafeno.

Pipeline:
1. Lê CIF com pymatgen.
2. Gera uma supercélula periódica 3x3x1.
3. Constrói grafo por distância.
4. Detecta ciclos/anéis.
5. Filtra anéis cujo centróide cai na célula central.
6. Projeta cada anel em 2D.
7. Aplica transformada de Joukowsky:
       w = z + a^2 / z
8. Calcula descritores:
       tamanho do anel,
       composição química,
       área,
       perímetro,
       anisotropia,
       área transformada,
       perímetro transformado,
       razão de área,
       razão de perímetro,
       anisotropia transformada.

Autor: exemplo para pesquisa / prototipagem.
"""

import argparse
import json
import math
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import networkx as nx

from pymatgen.core import Structure, Lattice
from pymatgen.core.periodic_table import Element


# ============================================================
# Utilidades geométricas
# ============================================================

def polygon_area_2d(points):
    """
    Área de um polígono 2D usando fórmula do sapateiro.
    points: array N x 2
    """
    points = np.asarray(points)
    x = points[:, 0]
    y = points[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def polygon_perimeter_2d(points):
    """
    Perímetro de um polígono 2D.
    """
    points = np.asarray(points)
    diffs = np.roll(points, -1, axis=0) - points
    return float(np.sum(np.linalg.norm(diffs, axis=1)))


def anisotropy_2d(points, eps=1e-12):
    """
    Anisotropia geométrica via razão dos autovalores da matriz de covariância.
    Para um anel quase circular/regular, tende a ~1.
    Para um anel alongado, cresce.
    """
    points = np.asarray(points)
    centered = points - points.mean(axis=0)
    cov = centered.T @ centered / max(len(points), 1)
    vals = np.linalg.eigvalsh(cov)
    vals = np.sort(vals)
    return float((vals[-1] + eps) / (vals[0] + eps))


def project_ring_to_2d(coords_3d):
    """
    Projeta coordenadas 3D do anel para o plano 2D principal via PCA/SVD.

    Isso torna o método aplicável mesmo quando a rede não está perfeitamente no plano xy.
    """
    coords_3d = np.asarray(coords_3d, dtype=float)
    center = coords_3d.mean(axis=0)
    X = coords_3d - center

    # SVD: os dois primeiros vetores definem o plano principal.
    _, _, vh = np.linalg.svd(X, full_matrices=False)
    basis = vh[:2].T

    coords_2d = X @ basis
    return coords_2d


def joukowsky_transform(points_2d, a=0.25, eps=1e-10):
    """
    Aplica transformada de Joukowsky aos pontos 2D do anel.

    z = x + i y
    w = z + a^2 / z

    Antes da transformação:
    - centraliza no centróide;
    - normaliza pelo raio médio.

    Retorna pontos transformados em 2D.
    """
    points_2d = np.asarray(points_2d, dtype=float)

    centered = points_2d - points_2d.mean(axis=0)
    z = centered[:, 0] + 1j * centered[:, 1]

    r = np.abs(z)
    r_mean = np.mean(r)

    if r_mean < eps:
        raise ValueError("Anel degenerado: raio médio muito pequeno.")

    z = z / r_mean

    # Evita singularidade em z = 0.
    z = np.where(np.abs(z) < eps, eps + 0j, z)

    w = z + (a ** 2) / z

    transformed = np.column_stack([w.real, w.imag])
    return transformed


# ============================================================
# Construção do grafo periódico
# ============================================================

def covalent_cutoff(site_i, site_j, scale=1.25, default_radius=0.75):
    """
    Cutoff químico baseado em raios covalentes.
    """
    el_i = Element(str(site_i.specie))
    el_j = Element(str(site_j.specie))

    r_i = el_i.average_cationic_radius or el_i.atomic_radius or default_radius
    r_j = el_j.average_cationic_radius or el_j.atomic_radius or default_radius

    return scale * float(r_i + r_j)


def build_periodic_2d_graph(
    structure,
    super_range=(-1, 0, 1),
    z_images=(0,),
    fixed_cutoff=None,
    cutoff_scale=1.25,
    max_bond=2.2,
):
    """
    Constrói grafo finito a partir de imagens periódicas.

    Para grafeno:
    - fixed_cutoff = 1.7 costuma funcionar bem.
    - super_range=(-1,0,1), z_images=(0,) gera 3x3x1.

    Cada nó é uma tupla:
        (índice_do_átomo_original, tx, ty, tz)

    onde tx, ty, tz são translações de célula.
    """

    lattice_matrix = structure.lattice.matrix
    nodes = []
    node_positions = {}
    node_species = {}
    node_frac = {}

    translations = []
    for tx in super_range:
        for ty in super_range:
            for tz in z_images:
                translations.append((tx, ty, tz))

    for tx, ty, tz in translations:
        shift_frac = np.array([tx, ty, tz], dtype=float)

        for i, site in enumerate(structure):
            node = (i, tx, ty, tz)

            frac = np.array(site.frac_coords, dtype=float) + shift_frac
            cart = frac @ lattice_matrix

            nodes.append(node)
            node_positions[node] = cart
            node_species[node] = str(site.specie)
            node_frac[node] = frac

    G = nx.Graph()
    for node in nodes:
        G.add_node(
            node,
            position=node_positions[node],
            species=node_species[node],
            frac=node_frac[node],
        )

    # Busca O(N^2). Para células pequenas/médias é suficiente.
    # Para sistemas grandes, trocar por KDTree.
    for idx_a in range(len(nodes)):
        node_a = nodes[idx_a]
        site_a = structure[node_a[0]]
        pos_a = node_positions[node_a]

        for idx_b in range(idx_a + 1, len(nodes)):
            node_b = nodes[idx_b]

            # Não conecta um átomo com sua própria imagem.
            if node_a[0] == node_b[0]:
                continue

            site_b = structure[node_b[0]]
            pos_b = node_positions[node_b]

            dist = np.linalg.norm(pos_a - pos_b)

            if fixed_cutoff is not None:
                cutoff = fixed_cutoff
            else:
                cutoff = covalent_cutoff(site_a, site_b, scale=cutoff_scale)

            cutoff = min(cutoff, max_bond)

            if 0.1 < dist <= cutoff:
                G.add_edge(node_a, node_b, distance=float(dist))

    return G


# ============================================================
# Detecção, ordenação e filtragem dos anéis
# ============================================================

def order_cycle_nodes(G, cycle):
    """
    Ordena os nós de um ciclo.

    networkx.minimum_cycle_basis retorna uma lista de nós do ciclo,
    mas nem sempre na ordem geométrica. Aqui reconstruímos a ordem
    usando o subgrafo induzido.

    Retorna None se o ciclo não for simples.
    """
    sub = G.subgraph(cycle)

    degrees = dict(sub.degree())
    if any(deg != 2 for deg in degrees.values()):
        return None

    start = cycle[0]
    ordered = [start]

    neighbors = list(sub.neighbors(start))
    prev = None
    current = start
    nxt = neighbors[0]

    while True:
        ordered.append(nxt)
        prev, current = current, nxt

        next_candidates = [n for n in sub.neighbors(current) if n != prev]

        if len(next_candidates) != 1:
            return None

        nxt = next_candidates[0]

        if nxt == start:
            break

        if nxt in ordered:
            return None

    return ordered


def canonical_cycle_key(cycle):
    """
    Chave canônica para remover duplicatas.

    Usa os nós completos, incluindo a imagem periódica.
    """
    return tuple(sorted(cycle))


def ring_centroid_fractional(G, cycle):
    """
    Centróide em coordenadas fracionárias expandidas.
    """
    fracs = np.array([G.nodes[n]["frac"] for n in cycle], dtype=float)
    return fracs.mean(axis=0)


def is_centroid_in_central_cell(frac_centroid, tol=1e-6):
    """
    Mantém anéis cujo centróide cai na célula original:
        0 <= x < 1
        0 <= y < 1

    Para sistemas 2D com vácuo em z, não impomos z.
    """
    x, y, _ = frac_centroid
    return (-tol <= x < 1.0 - tol) and (-tol <= y < 1.0 - tol)


def find_rings(
    G,
    min_size=3,
    max_size=12,
    central_only=True,
):
    """
    Detecta anéis usando minimum_cycle_basis.

    Para grafeno e redes 2D com anéis pequenos, isso funciona bem.
    Para redes muito complexas, pode-se substituir por algoritmos SSSR
    ou face-finding em grafos planares.
    """
    raw_cycles = nx.minimum_cycle_basis(G, weight="distance")

    rings = []
    seen = set()

    for cyc in raw_cycles:
        if not (min_size <= len(cyc) <= max_size):
            continue

        ordered = order_cycle_nodes(G, cyc)
        if ordered is None:
            continue

        key = canonical_cycle_key(ordered)
        if key in seen:
            continue
        seen.add(key)

        centroid_frac = ring_centroid_fractional(G, ordered)

        if central_only and not is_centroid_in_central_cell(centroid_frac):
            continue

        rings.append(ordered)

    return rings


# ============================================================
# Assinatura química do anel
# ============================================================

def minimal_rotation(seq):
    """
    Retorna rotação lexicograficamente mínima de uma sequência.
    """
    seq = list(seq)
    rotations = [tuple(seq[i:] + seq[:i]) for i in range(len(seq))]
    return min(rotations)


def canonical_species_signature(species):
    """
    Assinatura química do anel, independente do ponto inicial e do sentido.

    Exemplo:
        C-C-B-N-C-C
    """
    forward = minimal_rotation(species)
    backward = minimal_rotation(list(reversed(species)))
    best = min(forward, backward)
    return "-".join(best)


# ============================================================
# Descritor para cada anel
# ============================================================

def describe_ring(G, ring, joukowsky_a=0.25):
    """
    Calcula descritores geométricos e Joukowsky para um anel.
    """
    coords_3d = np.array([G.nodes[n]["position"] for n in ring], dtype=float)
    species = [G.nodes[n]["species"] for n in ring]

    coords_2d = project_ring_to_2d(coords_3d)

    area = polygon_area_2d(coords_2d)
    perimeter = polygon_perimeter_2d(coords_2d)
    aniso = anisotropy_2d(coords_2d)

    jw = joukowsky_transform(coords_2d, a=joukowsky_a)

    j_area = polygon_area_2d(jw)
    j_perimeter = polygon_perimeter_2d(jw)
    j_aniso = anisotropy_2d(jw)

    eps = 1e-12

    return {
        "ring_size": len(ring),
        "species_signature": canonical_species_signature(species),
        "species_count": dict(Counter(species)),

        "area": float(area),
        "perimeter": float(perimeter),
        "anisotropy": float(aniso),

        "j_area": float(j_area),
        "j_perimeter": float(j_perimeter),
        "j_anisotropy": float(j_aniso),

        "j_area_ratio": float(j_area / (area + eps)),
        "j_perimeter_ratio": float(j_perimeter / (perimeter + eps)),
    }


# ============================================================
# Agregação do descritor global
# ============================================================

def aggregate_descriptors(ring_descriptors):
    """
    Agrega descritores por tamanho de anel e por assinatura química.
    """
    if len(ring_descriptors) == 0:
        return {
            "n_rings": 0,
            "ring_size_counts": {},
            "ring_size_fractions": {},
            "by_ring_size": {},
            "by_species_signature": {},
        }

    df = pd.DataFrame(ring_descriptors)

    n_total = len(df)

    size_counts = df["ring_size"].value_counts().sort_index().to_dict()
    size_fractions = {
        int(k): float(v / n_total)
        for k, v in size_counts.items()
    }

    numeric_cols = [
        "area",
        "perimeter",
        "anisotropy",
        "j_area",
        "j_perimeter",
        "j_anisotropy",
        "j_area_ratio",
        "j_perimeter_ratio",
    ]

    by_size = {}
    for size, group in df.groupby("ring_size"):
        stats = {}
        for col in numeric_cols:
            stats[col + "_mean"] = float(group[col].mean())
            stats[col + "_std"] = float(group[col].std(ddof=0))
            stats[col + "_min"] = float(group[col].min())
            stats[col + "_max"] = float(group[col].max())
        stats["count"] = int(len(group))
        stats["fraction"] = float(len(group) / n_total)
        by_size[int(size)] = stats

    by_signature = {}
    for sig, group in df.groupby("species_signature"):
        stats = {}
        for col in numeric_cols:
            stats[col + "_mean"] = float(group[col].mean())
            stats[col + "_std"] = float(group[col].std(ddof=0))
        stats["count"] = int(len(group))
        stats["fraction"] = float(len(group) / n_total)
        by_signature[sig] = stats

    return {
        "n_rings": int(n_total),
        "ring_size_counts": {int(k): int(v) for k, v in size_counts.items()},
        "ring_size_fractions": size_fractions,
        "by_ring_size": by_size,
        "by_species_signature": by_signature,
    }


def compute_joukowsky_ring_descriptor(
    cif_file,
    fixed_cutoff=1.7,
    cutoff_scale=1.25,
    max_bond=2.2,
    min_ring_size=3,
    max_ring_size=12,
    joukowsky_a=0.25,
):
    """
    Função principal.
    """
    structure = Structure.from_file(cif_file)

    G = build_periodic_2d_graph(
        structure=structure,
        super_range=(-1, 0, 1),
        z_images=(0,),
        fixed_cutoff=fixed_cutoff,
        cutoff_scale=cutoff_scale,
        max_bond=max_bond,
    )

    rings = find_rings(
        G,
        min_size=min_ring_size,
        max_size=max_ring_size,
        central_only=True,
    )

    ring_descriptors = []
    for ring in rings:
        try:
            desc = describe_ring(G, ring, joukowsky_a=joukowsky_a)
            ring_descriptors.append(desc)
        except Exception as exc:
            print(f"[AVISO] Anel ignorado por erro geométrico: {exc}")

    global_descriptor = aggregate_descriptors(ring_descriptors)

    result = {
        "input_cif": cif_file,
        "n_atoms_unit_cell": len(structure),
        "formula": structure.composition.reduced_formula,
        "graph": {
            "n_nodes_supercell": G.number_of_nodes(),
            "n_edges_supercell": G.number_of_edges(),
        },
        "parameters": {
            "fixed_cutoff": fixed_cutoff,
            "cutoff_scale": cutoff_scale,
            "max_bond": max_bond,
            "min_ring_size": min_ring_size,
            "max_ring_size": max_ring_size,
            "joukowsky_a": joukowsky_a,
        },
        "descriptor": global_descriptor,
        "rings": ring_descriptors,
    }

    return result


# ============================================================
# CIF simples de grafeno para teste
# ============================================================

def write_graphene_cif(filename="graphene.cif", a=2.46, c=20.0):
    """
    Escreve um CIF simples para grafeno primitivo.

    Célula hexagonal:
        a = b = 2.46 Å
        gamma = 120°
        c grande para simular camada 2D.

    Dois átomos C por célula.
    """
    lattice = Lattice.hexagonal(a, c)

    species = ["C", "C"]

    frac_coords = [
        [1 / 3, 2 / 3, 0.5],
        [2 / 3, 1 / 3, 0.5],
    ]

    structure = Structure(
        lattice=lattice,
        species=species,
        coords=frac_coords,
        coords_are_cartesian=False,
    )

    structure.to(filename=filename)
    print(f"CIF de grafeno salvo em: {filename}")


# ============================================================
# CLI
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="Calcula descritor de anéis com transformada de Joukowsky para CIFs 2D."
    )

    parser.add_argument(
        "cif",
        nargs="?",
        help="Arquivo CIF de entrada."
    )

    parser.add_argument(
        "--make-graphene",
        metavar="OUT.cif",
        help="Gera um CIF simples de grafeno e sai."
    )

    parser.add_argument(
        "--cutoff",
        type=float,
        default=1.7,
        help="Cutoff fixo de ligação em Å. Para grafeno, 1.7 é adequado."
    )

    parser.add_argument(
        "--max-ring",
        type=int,
        default=12,
        help="Tamanho máximo de anel considerado."
    )

    parser.add_argument(
        "--joukowsky-a",
        type=float,
        default=0.25,
        help="Parâmetro a da transformada de Joukowsky."
    )

    parser.add_argument(
        "--json-out",
        default="descriptor.json",
        help="Arquivo JSON de saída."
    )

    parser.add_argument(
        "--csv-rings",
        default="rings.csv",
        help="Arquivo CSV com descritores individuais dos anéis."
    )

    args = parser.parse_args()

    if args.make_graphene:
        write_graphene_cif(args.make_graphene)
        return

    if args.cif is None:
        raise SystemExit("Forneça um arquivo CIF ou use --make-graphene.")

    result = compute_joukowsky_ring_descriptor(
        cif_file=args.cif,
        fixed_cutoff=args.cutoff,
        max_ring_size=args.max_ring,
        joukowsky_a=args.joukowsky_a,
    )

    with open(args.json_out, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, ensure_ascii=False)

    rings_df = pd.DataFrame(result["rings"])
    rings_df.to_csv(args.csv_rings, index=False)

    print("\n=== RESUMO ===")
    print(f"Arquivo: {result['input_cif']}")
    print(f"Fórmula: {result['formula']}")
    print(f"Átomos na célula unitária: {result['n_atoms_unit_cell']}")
    print(f"Nós no grafo 3x3x1: {result['graph']['n_nodes_supercell']}")
    print(f"Arestas no grafo 3x3x1: {result['graph']['n_edges_supercell']}")
    print(f"Número de anéis detectados na célula central: {result['descriptor']['n_rings']}")

    print("\nContagem por tamanho de anel:")
    for size, count in result["descriptor"]["ring_size_counts"].items():
        frac = result["descriptor"]["ring_size_fractions"][size]
        print(f"  anel {size}: {count}  fração={frac:.4f}")

    print("\nDescritores médios por tamanho de anel:")
    for size, stats in result["descriptor"]["by_ring_size"].items():
        print(f"\n  Ring size {size}")
        print(f"    count: {stats['count']}")
        print(f"    area_mean: {stats['area_mean']:.6f}")
        print(f"    perimeter_mean: {stats['perimeter_mean']:.6f}")
        print(f"    anisotropy_mean: {stats['anisotropy_mean']:.6f}")
        print(f"    j_area_ratio_mean: {stats['j_area_ratio_mean']:.6f}")
        print(f"    j_perimeter_ratio_mean: {stats['j_perimeter_ratio_mean']:.6f}")
        print(f"    j_anisotropy_mean: {stats['j_anisotropy_mean']:.6f}")

    print(f"\nJSON salvo em: {args.json_out}")
    print(f"CSV dos anéis salvo em: {args.csv_rings}")


if __name__ == "__main__":
    main()