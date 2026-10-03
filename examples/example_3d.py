from pathlib import Path

from periodic_ring_descriptors import compute_descriptor_3d, structure_feature_vector

structure = Path(__file__).parent / "structures" / "diamond.cif"
result = compute_descriptor_3d(
    structure,
    cutoff=1.895,
    max_ring_size=20,
    joukowsky_a=0.25,
    image_range=2,
)
features = structure_feature_vector(result)

print("3D descriptor")
print("atoms:", result["n_atoms"])
print("ring sizes:", result["ring_size_counts"])
print("rings per atom:", features["rings_per_atom"])
print("mean Joukowsky area ratio:", features["j_area_ratio_mean"])