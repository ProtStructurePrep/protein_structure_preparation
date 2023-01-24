.PHONY: install dev dev_environment test clean
NAME = hdmd

# check if we are in a conda environment
ifndef CONDA_PREFIX
$(error Please work in a conda environment)
endif

install:
	git checkout main; pip install -r requirements.txt; pip install .
dev:
ifeq (, $(shell which pre-commit))
	make dev_environment
endif
	git checkout develop
	pip uninstall $(NAME)
	pip install -e .
dev_environment:
	pip install -r requirements.txt
	pip install -r dev-requirements.txt
	pre-commit install; pre-commit autoupdate
test:
	pytest
jupyter:
	pip install ipykernel
	ipython kernel install --user --name=my-conda-env-kernel
clean:
	$(RM) -- *~ **/*~
	$(RM) -- *# **/*#
	$(RM) -r .cache .pytest_cache .libs tests/__pycache__ __pycache__ l3d.egg-info dist build $(NAME)/__pycache__
	pip uninstall $(NAME)
