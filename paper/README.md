# Reproducing the manuscript

This directory contains the exact source tables and scripts used to regenerate the figures for the PCCP manuscript.

## Requirements

Install the repository from its root with the paper extras:

```bash
python -m pip install -e ".[dev,paper]"
```

A TeX distribution providing `pdflatex` and `bibtex` is required. The optional `pdftotext` command enables the vector-text overlap audit.

## Reproduce everything

```bash
cd paper
bash reproduce.sh
```

The script:

1. regenerates all main and supplementary figures;
2. audits figure text for geometric overlap;
3. runs the descriptor regression/invariance tests from the repository root;
4. compiles the PCCP manuscript;
5. compiles the Electronic Supplementary Information.

Generated figure and PDF files are ignored by Git because they are deterministic build products.