.PHONY: test examples figures paper

test:
	python -m pytest -q

examples:
	python examples/example_2d.py
	python examples/example_3d.py

figures:
	cd paper && python scripts/make_figures.py

paper:
	cd paper && bash reproduce.sh