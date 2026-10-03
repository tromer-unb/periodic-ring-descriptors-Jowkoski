# Periodic Ring Descriptors for 2D and 3D Carbon Networks

[![tests](https://github.com/tromer-unb/periodic-ring-descriptors-Jowkoski/actions/workflows/tests.yml/badge.svg)](https://github.com/tromer-unb/periodic-ring-descriptors-Jowkoski/actions/workflows/tests.yml)
![Python](https://img.shields.io/badge/python-3.10%2B-blue)

This repository contains the reference implementation and paper-reproduction package for **periodic ring-geometry descriptors** designed to expose medium-range structural information in carbon networks.

The project provides two related descriptors:

- **2D descriptor:** periodic rings are represented by the bounded network motifs of planar carbon structures. The released implementation reproduces the paper calculations and was audited against explicit face tracing on all 120 allotropes used in the study.
- **3D descriptor:** rings are defined by a representative periodic bond plus every shortest alternate path between its endpoints. Translationally equivalent cycles are canonicalized, avoiding the crystallographic-cell dependence of a finite cycle basis.

For each ring, the code reports topology, projected geometry, and an optional orientation-averaged Joukowsky response.

> **Research status.** This is the reference research code accompanying the manuscript. The descriptor definition, validation tests, and paper-reproduction data are versioned together so that published numerical results can be traced to executable code.

## Installation

Clone the repository and install it in editable mode:

```bash
git clone https://github.com/tromer-unb/periodic-ring-descriptors-Jowkoski.git
cd periodic-ring-descriptors-Jowkoski
python -m pip install -e ".[dev,paper]"
```

The core package depends on NumPy, pandas, NetworkX, and pymatgen. Matplotlib is only required to reproduce the paper figures.

## Quick start

### 2D

```python
from periodic_ring_descriptors import compute_descriptor_2d, structure_feature_vector

result = compute_descriptor_2d(
    "structure.cif",
    cutoff=1.895,
    max_ring_size=12,
    a=0.25,
)
features = structure_feature_vector(result)

print(result["descriptor"]["ring_size_counts"])
print(features["area_mean"])
print(features["j_area_ratio_mean"])
```

### 3D

```python
from periodic_ring_descriptors import compute_descriptor_3d, structure_feature_vector

result = compute_descriptor_3d(
    "structure.cif",
    cutoff=1.895,
    max_ring_size=20,
    joukowsky_a=0.25,
    image_range=2,
)
features = structure_feature_vector(result)

print(result["ring_size_counts"])
print(result["rings_per_atom"])
print(features["j_area_ratio_mean"])
```

Runnable examples are provided in [`examples/`](examples/).

## Command-line interface

After installation:

```bash
ringdesc 2d examples/structures/graphene.cif --cutoff 1.895 --rings-out rings_2d.csv
ringdesc 3d examples/structures/diamond.cif --cutoff 1.895 --image-range 2 --rings-out rings_3d.csv
```

Each command writes a JSON descriptor. The optional CSV contains one row per detected ring.

## Descriptor content

Each ring is projected onto its local principal plane and characterized by:

- ring size and chemical composition;
- projected area `A`, perimeter `P`, mean radius, and anisotropy;
- orientation-averaged transformed area `A_J`, perimeter `P_J`, and anisotropy;
- dimensionless response ratios `R_A = A_J/A` and `R_P = P_J/P`.

The normalized, scale-preserving Joukowsky map is

`w = r [ u + a^2/u ]`, with `u = z/r`,

evaluated after aligning each ring edge with the positive real axis and averaging the resulting geometric measures over all edge alignments.

For a flat, model-ready dictionary, use `structure_feature_vector(result)`. The helper summarizes the common topology/geometry/Joukowsky blocks but intentionally does **not** add dataset-specific defect labels or chemistry metadata used in the paper's machine-learning benchmarks.

## Validation

The repository includes explicit regression and invariance tests covering:

- 2D physical units and the Joukowsky identity limit;
- rotation, reflection, cycle direction, and vertex-order invariance;
- graphene primitive-cell/supercell consistency;
- diamond primitive/conventional/supercell consistency;
- 3D origin shifts, reflections, and site reordering;
- lonsdaleite cell-replication consistency.

Run:

```bash
pytest -q
```

The GitHub Actions workflow runs the test suite on Python 3.10, 3.11, and 3.12.

## Reproducing the paper

The [`paper/`](paper/) directory contains the PCCP LaTeX source, Supplementary Information, aggregated data used in every figure, and the plotting/audit scripts.

```bash
cd paper
bash reproduce.sh
```

The script regenerates all main and supplementary figures, checks for text overlap in the vector figures when `pdftotext` is available, runs the descriptor tests, and compiles the manuscript and ESI.

The repository intentionally contains the **aggregated numerical data required for figure reproduction**, not thousands of raw relaxation files. This keeps the research artifact auditable and lightweight. See [`docs/reproducibility.md`](docs/reproducibility.md).

## Repository layout

```text
src/periodic_ring_descriptors/   installable descriptor package
tests/                           invariance and regression tests
examples/                        minimal 2D/3D usage examples
docs/                            scientific and API documentation
paper/                           manuscript + ESI + figure data/scripts
.github/workflows/               continuous integration
```

## Scientific scope and caveats

The bond graph is part of the descriptor definition. The paper uses a 1.895 Å cutoff for the reported carbon benchmarks and explicitly audits cutoff sensitivity. For new materials or chemistry, users should establish a physically defensible connectivity criterion rather than copying this value blindly.

The 3D descriptor uses an **edge-shortest periodic ring definition**. It is tested for the crystal representations used in the paper; the code does not claim that every possible ring convention is equivalent in arbitrary periodic networks.

The Joukowsky layer should be treated as an optional shape-response block. In the paper it is largely redundant with raw geometry for planar 2D carbon, but contributes additional information in the corrected 3D benchmark.

## Documentation

- [2D descriptor](docs/descriptor_2d.md)
- [3D descriptor](docs/descriptor_3d.md)
- [Output and feature schema](docs/output_schema.md)
- [Paper reproducibility](docs/reproducibility.md)
- [Scientific design notes](docs/scientific_design.md)

## Citation

The manuscript is currently being prepared for submission to **Physical Chemistry Chemical Physics (PCCP)**. A final citation and archival DOI will be added here when available.

Repository: <https://github.com/tromer-unb/periodic-ring-descriptors-Jowkoski>