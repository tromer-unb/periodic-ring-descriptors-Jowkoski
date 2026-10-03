#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
cd "$HERE"

echo "[1/5] Rebuilding figures"
python scripts/make_figures.py

echo "[2/5] Auditing figure text layout"
if command -v pdftotext >/dev/null 2>&1; then
  python scripts/audit_figure_text.py
else
  echo "pdftotext not available; skipping geometric text-overlap audit."
fi

echo "[3/5] Running descriptor tests"
(
  cd "$ROOT"
  python -m pytest -q
)

echo "[4/5] Compiling PCCP manuscript"
pdflatex -interaction=nonstopmode main.tex >/dev/null
bibtex main >/dev/null
pdflatex -interaction=nonstopmode main.tex >/dev/null
pdflatex -interaction=nonstopmode main.tex >/dev/null

echo "[5/5] Compiling Electronic Supplementary Information"
pdflatex -interaction=nonstopmode supplementary.tex >/dev/null
pdflatex -interaction=nonstopmode supplementary.tex >/dev/null

echo "Reproduction completed successfully."