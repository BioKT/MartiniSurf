<p align="center">
  <img src="logo_martini.png" alt="MartiniSurf Logo">
</p>

<h1 align="center">MartiniSurf</h1>

<p align="center">
A Python toolkit for automatic protein orientation, surface generation,  
and Gō–Martini coarse-grained simulation setup.
</p>


<p align="center">

  <!-- Build status -->
  <a href="https://github.com/jjimenezgar/MartiniSurf/actions/workflows/python-ci.yml">
    <img src="https://github.com/jjimenezgar/MartiniSurf/actions/workflows/python-ci.yml/badge.svg" alt="CI">
  </a>

  <!-- Code style -->

  <!-- flake8 -->
  <img src="https://img.shields.io/badge/lint-flake8-blue" alt="flake8">

  <!-- ruff -->
  <img src="https://img.shields.io/badge/lint-ruff-red" alt="ruff">

  <!-- mypy -->
  <img src="https://img.shields.io/badge/types-mypy-green" alt="mypy">

  <!-- Pytest -->
  <img src="https://img.shields.io/badge/tests-passing-brightgreen" alt="tests">

  <!-- Coverage -->
  <a href="https://codecov.io/gh/jjimenezgar/MartiniSurf">
    <img src="https://codecov.io/gh/jjimenezgar/MartiniSurf/branch/master/graph/badge.svg" alt="Coverage">
  </a>

  <!-- Python versions -->
  <img src="https://img.shields.io/badge/python-3.9%20|%203.10%20|%203.11-blue.svg">
  
  <!-- Google Colab -->
  <a href="https://colab.research.google.com/github/jjimenezgar/MartiniSurf/blob/master/martinisurf/examples/colab_demo.ipynb">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab">
  </a>

</p>


---

# 🧬 Overview

**MartiniSurf** is a complete, automated pipeline for preparing **Martini 3**,  
**Gō–Martini**, and **Martini2-DNA** simulation systems involving:

✔️ Proteins  
✔️ DNA (single- and double-stranded)  
✔️ Hybrid protein–surface or DNA–surface systems  

It provides a modular yet fully automatic way to:

- ⚙️ **Generate coarse-grained models** using  
  - `martinize2` (proteins)  
  - `martinize-dna.py` (DNA, Python-2 compatible)
- 🧱 **Build Martini surfaces** (custom bead type, charge, spacing)
- 📐 **Orient biomolecules** on surfaces via multi-anchor geometry
- 🔗 **Define unlimited anchor groups**
- 🔒 **Generate position restraints** for Gō–Martini anchoring
- 🧰 **Produce full GROMACS-ready directories**
  - Minimization, NVT, NPT, production `.mdp` files

---


## 🔧 **Installation**

- conda create -n martinisurf
- conda activate martinisurf
- pip install -r surfmartini
- conda install pip
- pip install .

## ⚡ **Quickstart**
### 🧬 Protein-Support
martinisurf \
    --pdb 1RJW \
    --moltype Protein \
    --lx 20 --ly 20 \
    --surface-bead C1 \
    --anchor 1 8 10 11 \
    --anchor 2 1025 1027 1028 \
    --dist 10

### 🔍 What this does

- 📥 Downloads **1RJW** from RCSB  
- 🧬 Runs **martinize2** to generate the coarse-grained model  
- 🧱 Builds a **20 × 20 nm** Martini surface  
- 📐 Orients the enzyme using **two anchor groups**  
- 🧰 Generates complete **GROMACS-ready simulation files**  


### 🧬 DNA-Support

martinisurf \
    --dna \
    --pdb 4C64 \
    --moltype DNA \
    --dnatype ds-stiff \
    --lx 5 --ly 5 \
    --surface-bead C1 \
    --anchor 1 1 \
    --anchor 2 24 \
    --dist 10

### 🔍 What this does

- 📥 Downloads **4C64** from RCSB  
- 🧬 Runs **martinize-dna.py** to generate the coarse-grained model  
- 🧱 Builds a **5 × 5 nm** Martini surface  
- 📐 Orients the enzyme using **two anchor groups**  
- 🧰 Generates complete **GROMACS-ready simulation files**  

---

📜 Third-Party Software and Licensing Notice

MartiniSurf relies on several external open-source tools and scientific libraries.
These components are not included in this repository (unless explicitly stated) and remain governed by their original licenses.

Users must ensure that they comply with the licenses of all dependencies when using MartiniSurf.

🔧 Core Python Dependencies

MartiniSurf depends on the following libraries:

Package	License	Purpose
NumPy	        BSD-3-Clause	Scientific computing
SciPy	        BSD-3-Clause	Numerical routines
MDTraj	        BSD-3-Clause	Trajectory handling
MDAnalysis	LGPL-2.1	Atomistic topology & trajectory analysis
PyVista	MIT	        3D visualization
Vermouth/martinize2	Apache-2.0	Martini CG model construction


🧬 Martini and DNA Martini Tools

MartiniSurf interfaces with external Martini tooling, including:

martinize2 (Apache-2.0), used for protein coarse-graining

martinize-dna.py, originally developed within the Marrink Lab, distributed under research-friendly open licensing (GPL-based in earlier versions)

⚠️ These tools are not distributed within MartiniSurf.
They are invoked only when present on the user’s system, and their licenses remain fully applicable.
