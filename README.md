[Installation](#installation) |  [Requirements](#requirements) | [Example](#example) | [Usage](#usage) | [Command line](#command-line) | [Acknowledgment](#acknowledgment)

Protein-Ligand Complex MD simulations
=================
This repository contains the basic steps for parsing proteins and selecting its ligands. Also, it contains the steps of preparing the protein and ligands, by fixing errors in the PDB file.

Table of contents
=================

* [Installation](#installation)
* [Requirements](#requirements)
* [Contributing](#contributing)
* [Acknowledgment](#acknowledgment)


Installation
============
To install ....

- Available on `conda-forge` channel

    ```bash
    conda install nglview -c conda-forge
    # might need: jupyter-nbextension enable nglview --py --sys-prefix

    # if you already installed nglview, you can `upgrade`
    conda upgrade nglview --force
    # might need: jupyter-nbextension enable nglview --py --sys-prefix
    ```

Requirements
=================
* [MDTraj](https://github.com/mdtraj/mdtraj)
* [nglview](https://github.com/nglviewer/nglview)
* [OpenMM](https://github.com/openmm/openmm) and [OpenMM Forcefields](https://github.com/openmm/openmmforcefields)
* [RDKit](https://github.com/rdkit/rdkit)
* [PDBFixer](https://github.com/openmm/pdbfixer).


