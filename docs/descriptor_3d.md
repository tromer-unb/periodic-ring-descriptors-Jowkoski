# 3D periodic edge-shortest ring descriptor

## Why a different ring definition is required

A finite graph cycle basis is not a crystallographic invariant. In the original 3D prototype, the same diamond lattice produced incompatible cycle sets depending on whether it was represented by a primitive cell or a larger supercell.

The released 3D descriptor therefore defines rings directly from periodic bonds rather than from a finite cycle basis.

## Edge-shortest periodic rings

For each representative periodic bond whose midpoint lies in the reference cell:

1. remove the bond temporarily;
2. find every shortest alternate path between its endpoints in the periodic image graph;
3. combine each shortest path with the removed bond;
4. canonicalize the cycle under starting vertex, traversal direction, and uniform lattice translation;
5. remove translational duplicates.

This construction is tested in the repository for diamond and lonsdaleite cell replication, origin shifts, rigid transformations, and site reordering.

## Recommended call

```python
from periodic_ring_descriptors import compute_descriptor_3d

result = compute_descriptor_3d(
    "diamond.cif",
    cutoff=1.895,
    max_ring_size=20,
    joukowsky_a=0.25,
    image_range=2,
)
```
## Key output

The 3D result includes:

- `ring_size_counts`
- `rings_per_atom`
- `ring_memberships_per_atom`
- `n_central_bonds`
- coordination range for central-cell atoms
- complete per-ring geometric and Joukowsky quantities

For ideal cubic diamond, the validated reference result contains only six-membered rings, with 2 rings per atom and 12 ring memberships per atom.

## `image_range`

`image_range` controls the explicit periodic image graph. A value of 1 is sufficient for the 64-atom defect supercells used in the paper. Smaller primitive cells used in invariance tests may require 2 so that every relevant alternate path is represented.

## Scope

The implementation establishes representation invariance for the tested carbon crystals; it does not assert that all ring definitions become equivalent in arbitrary periodic networks. If the bonding topology or expected ring sizes differ substantially from the paper systems, convergence with respect to cutoff, image range, and maximum ring size should be checked explicitly.