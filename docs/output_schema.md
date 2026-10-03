# Output and feature schema

Both descriptor functions return a Python dictionary and expose the full list of detected rings under `result["rings"]`.

## Per-ring fields

| Field | Meaning | Units |
|---|---|---|
| `ring_size` | number of vertices | dimensionless |
| `species_signature` | order-invariant chemical signature | — |
| `area` | projected polygon area | Å² |
| `perimeter` | projected polygon perimeter | Å |
| `anisotropy` | covariance eigenvalue ratio | dimensionless |
| `mean_radius` | mean radial distance in local plane | Å |
| `j_area` | orientation-averaged transformed area | Å² |
| `j_perimeter` | orientation-averaged transformed perimeter | Å |
| `j_anisotropy` | transformed anisotropy | dimensionless |
| `j_area_ratio` | `j_area / area` | dimensionless |
| `j_perimeter_ratio` | `j_perimeter / perimeter` | dimensionless |

## Flat model-ready features

`structure_feature_vector(result)` returns a flat dictionary containing ring-count statistics, ring-size fractions, and mean/standard deviation of the common raw and Joukowsky geometry fields.

The helper is intentionally chemistry-agnostic. The paper's defect benchmarks additionally use dataset-specific variables such as defect count and B/N ring composition; those benchmark tables are reproduced from the archived paper data rather than silently injected into the public descriptor API.