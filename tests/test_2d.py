from pathlib import Path

import numpy as np
import pytest

from periodic_ring_descriptors._core_2d import polygon_area_2d, polygon_perimeter_2d
from periodic_ring_descriptors import compute_descriptor_2d, joukowsky_metrics_v2

ROOT = Path(__file__).resolve().parent / "structures"


@pytest.fixture
def hexagon():
    angles = np.arange(6) * np.pi / 3
    return 1.42 * np.column_stack((np.cos(angles), np.sin(angles)))


def test_identity_has_physical_units(hexagon):
    result = joukowsky_metrics_v2(hexagon, a=0)
    assert np.isclose(result["j_area"], polygon_area_2d(hexagon))
    assert np.isclose(result["j_perimeter"], polygon_perimeter_2d(hexagon))


def test_rotation_reflection_and_vertex_order(hexagon):
    angle = 0.371
    rotation = np.array([[np.cos(angle), -np.sin(angle)],
                         [np.sin(angle), np.cos(angle)]])
    irregular = np.array([[0., 0.], [1.3, 0.1], [1.7, 1.2],
                          [0.4, 1.8], [-0.5, 0.9]])
    for shape in (hexagon, irregular):
        expected = joukowsky_metrics_v2(shape)
        variants = [shape @ rotation.T, shape[::-1], shape[:, ::-1], np.roll(shape, 2, axis=0)]
        for points in variants:
            actual = joukowsky_metrics_v2(points)
            for key in ("j_area", "j_perimeter", "j_anisotropy"):
                assert np.isclose(actual[key], expected[key], atol=1e-10)


def test_length_scaling_preserves_geometry_rules(hexagon):
    original = joukowsky_metrics_v2(hexagon)
    doubled = joukowsky_metrics_v2(2 * hexagon)
    assert np.isclose(doubled["j_area"], 4 * original["j_area"])
    assert np.isclose(doubled["j_perimeter"], 2 * original["j_perimeter"])
    assert np.isclose(doubled["j_anisotropy"], original["j_anisotropy"])


def test_vertex_at_center_is_rejected():
    points = np.array([[0., 0.], [1., 0.], [0., 1.], [-1., -1.]])
    with pytest.raises(ValueError, match="singular"):
        joukowsky_metrics_v2(points)


def test_graphene_reference_values():
    result = compute_descriptor_2d(ROOT / "graphene.cif")
    assert result["descriptor_version"] == 2
    assert result["descriptor"]["ring_size_counts"] == {6: 1}
    ring = result["rings"][0]
    assert np.isclose(ring["j_area_ratio"], 0.99609375, atol=1e-7)
    assert np.isclose(ring["j_perimeter_ratio"], 1.001007279, atol=1e-7)


def test_graphene_cell_replication_invariance():
    primitive = compute_descriptor_2d(ROOT / "graphene.cif")
    supercell = compute_descriptor_2d(ROOT / "graphene_supercell.cif")
    assert primitive["descriptor"]["ring_size_counts"] == {6: 1}
    assert supercell["descriptor"]["ring_size_counts"] == {6: 16}
    assert np.isclose(
        len(primitive["rings"]) / primitive["n_atoms_unit_cell"],
        len(supercell["rings"]) / supercell["n_atoms_unit_cell"],
    )
    for key in (
        "area", "perimeter", "anisotropy",
        "j_area", "j_perimeter", "j_anisotropy",
        "j_area_ratio", "j_perimeter_ratio",
    ):
        left = np.mean([ring[key] for ring in primitive["rings"]])
        right = np.mean([ring[key] for ring in supercell["rings"]])
        assert np.isclose(left, right, atol=1e-8), key