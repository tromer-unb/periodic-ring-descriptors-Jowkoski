from pathlib import Path

import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from periodic_ring_descriptors import compute_descriptor_3d_structure

ROOT = Path(__file__).resolve().parent / "structures"


def _diamond_cells():
    source = Structure.from_file(ROOT / "diamond_reference.cif")
    analyzer = SpacegroupAnalyzer(source)
    primitive = analyzer.get_primitive_standard_structure()
    conventional = analyzer.get_conventional_standard_structure()
    supercell = primitive.copy()
    supercell.make_supercell([2, 2, 2])
    return primitive, conventional, supercell


def _summary(structure, image_range):
    result = compute_descriptor_3d_structure(
        structure, cutoff=1.895, max_ring_size=12,
        joukowsky_a=0.25, image_range=image_range,
    )
    ratios = np.array([ring["j_area_ratio"] for ring in result["rings"]])
    return result, float(ratios.mean())


def _fingerprint(structure, image_range=2):
    result = compute_descriptor_3d_structure(
        structure, cutoff=1.895, max_ring_size=20,
        joukowsky_a=0.25, image_range=image_range,
    )
    keys = (
        "area", "perimeter", "anisotropy", "mean_radius",
        "j_area", "j_perimeter", "j_anisotropy",
        "j_area_ratio", "j_perimeter_ratio",
    )
    fp = {
        "ring_size_counts_per_atom": {
            k: v / len(structure) for k, v in result["ring_size_counts"].items()
        },
        "rings_per_atom": result["rings_per_atom"],
        "ring_memberships_per_atom": result["ring_memberships_per_atom"],
    }
    for key in keys:
        fp[key] = float(np.mean([ring[key] for ring in result["rings"]]))
    return fp


def _assert_fp_equal(left, right):
    assert left["ring_size_counts_per_atom"] == right["ring_size_counts_per_atom"]
    for key in left:
        if key != "ring_size_counts_per_atom":
            assert np.isclose(left[key], right[key], atol=1e-10), key


def test_diamond_cell_invariance():
    primitive, conventional, supercell = _diamond_cells()
    cases = [_summary(primitive, 2), _summary(conventional, 2), _summary(supercell, 1)]
    reference_ratio = cases[0][1]
    for result, ratio in cases:
        assert result["ring_size_counts"] == {6: result["n_rings"]}
        assert np.isclose(result["rings_per_atom"], 2.0)
        assert np.isclose(result["ring_memberships_per_atom"], 12.0)
        assert np.isclose(ratio, reference_ratio, atol=1e-12)


def test_diamond_coordination():
    primitive, _, _ = _diamond_cells()
    result, _ = _summary(primitive, 2)
    assert result["degree_min"] == 4
    assert result["degree_max"] == 4
    assert result["n_central_bonds"] == 4


def test_origin_shift_rotation_reflection_and_site_order():
    _, conventional, _ = _diamond_cells()
    ref = _fingerprint(conventional, 2)

    shifted = Structure(
        conventional.lattice, conventional.species,
        (conventional.frac_coords + np.array([0.173, 0.287, 0.419])) % 1.0,
        coords_are_cartesian=False,
    )
    _assert_fp_equal(ref, _fingerprint(shifted, 2))

    angle = 0.731
    axis = np.array([1.0, 2.0, -0.5])
    axis /= np.linalg.norm(axis)
    x, y, z = axis
    c, sn = np.cos(angle), np.sin(angle)
    C = 1.0 - c
    rotation = np.array([
        [c+x*x*C, x*y*C-z*sn, x*z*C+y*sn],
        [y*x*C+z*sn, c+y*y*C, y*z*C-x*sn],
        [z*x*C-y*sn, z*y*C+x*sn, c+z*z*C],
    ])
    rotated = Structure(
        Lattice(conventional.lattice.matrix @ rotation.T),
        conventional.species, conventional.frac_coords,
        coords_are_cartesian=False,
    )
    _assert_fp_equal(ref, _fingerprint(rotated, 2))

    reflection = np.diag([-1.0, 1.0, 1.0])
    mirrored = Structure(
        Lattice(conventional.lattice.matrix @ reflection.T),
        conventional.species, conventional.frac_coords,
        coords_are_cartesian=False,
    )
    _assert_fp_equal(ref, _fingerprint(mirrored, 2))

    order = list(reversed(range(len(conventional))))
    reordered = Structure(
        conventional.lattice,
        [conventional[i].specie for i in order],
        [conventional[i].frac_coords for i in order],
        coords_are_cartesian=False,
    )
    _assert_fp_equal(ref, _fingerprint(reordered, 2))


def test_lonsdaleite_cell_invariance():
    lonsdaleite = Structure.from_file(ROOT / "lonsdaleite.cif")
    supercell = lonsdaleite.copy()
    supercell.make_supercell([2, 2, 2])
    ref = _fingerprint(lonsdaleite, 2)
    sc = _fingerprint(supercell, 1)
    _assert_fp_equal(ref, sc)
    assert ref["ring_size_counts_per_atom"] == {6: 2.0}
    assert np.isclose(ref["rings_per_atom"], 2.0)