###
### Task runner for code development
###

.PHONY: test

test:
	@make pytest
	@make doctest
	@make lint
	@make format-check

.PHONY: pytest
pytest:
	python -m pytest

.PHONY: doctest
doctest:
	sphinx-build -b doctest docs/source docs/build/doctest

.PHONY: lint
lint:
	ruff check .

.PHONY: format
format:
	ruff format .

.PHONY: format-check
format-check:
	ruff format --check .

.PHONY: docs
docs:
	sphinx-build -b html docs/source docs/build/html

.PHONY: clean
clean:
	rm -rf docs/build .pytest_cache .ruff_cache build
	find . -type d -name "*.egg-info" -exec rm -rf {} +
	find . -type d -name "__pycache__" -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete
