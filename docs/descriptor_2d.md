# 2D periodic ring descriptor

## Definition

The 2D descriptor operates on a periodic bond graph generated from a planar structure. For the published carbon benchmark, the conceptual ring definition is the bounded face of the periodic planar graph.

The released reference implementation reproduces the historical paper calculations using a distance-weighted minimum cycle basis with central-cell filtering. This implementation was audited against explicit planar face tracing on all 120 2D allotropes used in the study: ring counts, complete ring-size multisets, and structure-level geometric means were identical to numerical precision.

## Ring geometry

For each ring, the atomic coordinates are unwrapped across periodic boundaries and projected onto the first two principal directions. The code reports:

- `ring_size`
- `species_signature` and `species_count`
- `area`
- `perimeter`
- `anisotropy`
- `mean_radius`
- `j_area`, `j_perimeter`, `j_anisotropy`
- `j_area_ratio`, `j_perimeter_ratio`

Area is reported in Å² and perimeter/radius in Å when the input structure uses Ångström coordinates. The two Joukowsky ratios are dimensionless.
## Recommended call

```python
from periodic_ring_descriptors import compute_descriptor_2d

result = compute_descriptor_2d(
    "graphene.cif",
    cutoff=1.895,
    max_ring_size=12,
    a=0.25,
)
```

### Parameters

- `cutoff`: fixed bond cutoff in Å. The paper uses 1.895 Å.
- `max_ring_size`: maximum accepted ring size. The 2D paper analysis uses 12 for the main descriptor extraction; larger rings can be retained if required by a new dataset.
- `a`: dimensionless Joukowsky parameter. The paper uses 0.25 and reports a sensitivity study.

## Important scope note

The 1.895 Å cutoff is justified for the published 120-allotrope dataset by a global neighbor gap: the largest third-neighbor distance is below the smallest fourth-neighbor distance. That observation is specific to the reported carbon dataset. For a new system, bond connectivity should be audited independently.