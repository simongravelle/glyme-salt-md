[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20411031.svg)](https://doi.org/10.5281/zenodo.20411031)

# Cs-TFSI-PEO molecular system builder

This repository contains tools to generate atomistic molecular dynamics simulation systems of polyethylene oxide (PEO) electrolytes with CsTFSI salt.

The code builds polymer chains, inserts salt ions into a simulation box, and outputs GROMACS-compatible topology and coordinate files for bulk simulations.

The scripts in this repository are associated with the study: “Molecular-to-polymeric crossover in ion diffusion in glyme-based electrolytes: from vehicular to hopping transport” by Aicha Jani, Simon Gravelle, Mehdi Zeghal, Pawel Wzietek, and Patrick Judeinstein.

Input configurations for different polymer chain lengths are provided in the repository, for n = 1 to 88 (`PEO_1_6000`, `PEO_2_6000`, `PEO_4_6000`, `PEO_8_6000`, `PEO_24_6000`, `PEO_44_6000`, `PEO_88_6000`).

## Repository structure

The most important codes to regenerate the input files and run the GROMACS simulation are contained in the `create-system/` folder:

- `create_salt_peo_boxes.py`  
  Main script to generate bulk PEO–CsTFSI systems and write GROMACS input files (`.gro`, `.top`).

- `peo_generator.py`  
  Generates single PEO chains and their bonded topology (atoms, bonds, angles, dihedrals).

- `utilities.py`  
  Helper functions for random insertion and overlap detection.

- `DATA/`  
  Contains structural building blocks and force-field input data:
  - PEO monomer and end groups
  - TFSI molecular structure and bonded topology
  - CHARMM-related parameter files

- `ff/`  
  Force field definitions for:
  - PEO polymer
  - TFSI anion
  - Cs⁺ cation
  - Generic force-field parameters

- `input/`  
  Example GROMACS `.mdp` files for:
  - Energy minimization
  - NVT equilibration
  - NPT equilibration
  - Production runs
