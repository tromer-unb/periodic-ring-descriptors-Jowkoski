# Paper reproducibility

The `paper/` directory is a self-contained research artifact for the manuscript.

## Included

- PCCP LaTeX source for the main paper and ESI;
- RSC bibliography style and BibTeX database;
- aggregated CSV/JSON values underlying every figure;
- figure-generation script;
- vector-text overlap audit;
- descriptor regression/invariance tests through the repository root.

## Reproduce

From a configured checkout:

```bash
python -m pip install -e ".[dev,paper]"
cd paper
bash reproduce.sh
```

A standard TeX installation with `pdflatex` and `bibtex` is required to compile the manuscript. `pdftotext` is optional and is used only for the automated vector-text overlap audit.

## What is intentionally not duplicated

The repository does not include thousands of raw MACE-relaxed structures or large intermediate model objects. The paper figures are reproduced from the exact aggregated data used in the manuscript. Raw structural datasets and relaxation records can be archived separately when the manuscript is submitted/published.

## Provenance

The descriptor source, tests, paper source, and figure tables are versioned in the same Git commit. This makes a paper revision traceable to the exact code and numerical tables used to render it.