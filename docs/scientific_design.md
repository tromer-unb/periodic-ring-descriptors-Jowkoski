# Scientific design notes

## Three information layers

The descriptor separates three feature families:

1. **Topology** — ring sizes, counts, fractions, and coordination/network statistics.
2. **Raw ring geometry** — projected area, perimeter, radius, and anisotropy.
3. **Joukowsky response** — transformed geometry and dimensionless response ratios.

This separation is deliberate. It allows the transformed layer to be tested by ablation rather than treated as an inseparable black-box feature generator.

## Joukowsky response

After centering a projected ring, the complex coordinate is normalized by its mean radius. For every edge orientation, the ring is rotated so that the edge points along the positive real axis, the scale-preserving map `w = r (u + a²/u)` is applied, and transformed geometric measures are averaged over edge alignments.

The orientation averaging removes dependence on the arbitrary in-plane axis chosen by the local projection.

## What the paper finds

- In planar 2D carbon, raw and transformed geometry are highly correlated; Joukowsky features do not provide a robust incremental gain once raw geometry is present.
- In the corrected 3D benchmark, the transformed block adds measurable information to the raw ring model.
- Local descriptors and ring descriptors dominate in different structural regimes, and hybrid models show that the two scales can be complementary.

These are empirical findings for the reported datasets, not universal guarantees.