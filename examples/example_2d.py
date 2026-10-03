from pathlib import Path

from periodic_ring_descriptors import compute_descriptor_2d, structure_feature_vector

structure = Path(__file__).parent / "structures" / "graphene.cif"
result = compute_descriptor_2d(
    structure,
    cutoff=1.895,
    max_ring_size=12,
    a=0.25,
)
features = structure_feature_vector(result)

print("2D descriptor")
print("atoms:", result["n_atoms_unit_cell"])
print("ring sizes:", result["descriptor"]["ring_size_counts"])
print("mean area:", features["area_mean"])
print("mean Joukowsky area ratio:", features["j_area_ratio_mean"])