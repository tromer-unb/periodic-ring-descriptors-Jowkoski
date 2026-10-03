# Contributing

This repository is maintained as reference research software. Contributions that improve correctness, reproducibility, documentation, or computational efficiency are welcome.

## Before opening a pull request

1. Create a focused branch.
2. Install the development environment with `python -m pip install -e ".[dev,paper]"`.
3. Run `pytest -q`.
4. If figures change, run `cd paper && python scripts/make_figures.py` and the text-overlap audit.
5. Describe whether the change affects the descriptor definition, only the implementation, or only documentation.

Changes to the ring definition, periodic canonicalization, default cutoff, Joukowsky normalization, or aggregation semantics should include a regression test and a short scientific rationale.

Please do not mix scientific-definition changes with large formatting refactors in the same pull request.