<p align="center">
  <img src="logo_martini.png" alt="MartiniSurf Logo" width="240">
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
</p>




---

# 🧬 Overview

**MartiniSurf** is a high-level toolkit that automates the complete workflow
required to simulate enzyme immobilization using **Martini 3** and **Gō–Martini**.

It provides:

- ⚙️ **Automated construction of hexagonal surfaces** (GRO + ITP)
- 📐 **Automatic protein orientation** on surfaces using anchors
- 🔗 **Multi-anchor support** (Custom residues)
- 🧱 **Integration with martinize3 and Gō–Martini**
- 🧰 **Full pipeline mode:**  
  Generates *all* topology, index, restraints, and MDP files
- 🧬 **Python API** + **command-line interface**
- 🚀 Reproducible simulation folders for HPC clusters (SLURM-ready)


---

# ✨ Features

- 🧱 Build Martini-compatible **hex-lattice surfaces**
- 🧬 **Convert and orient enzymes** using PDB/GRO inputs
- ⛓️ Define **anchor residues** for directional immobilization
- 📦 Create **GROMACS-ready** directories:
  - `0_topology/`
  - `1_mdp/`
  - `2_system/`
- 🔧 Generate:
  - `system.top`
  - `system_res.top`
  - `Active.itp`
  - `Active_res.itp`
  - `surface.itp`
  - `index.ndx`
  - all `.mdp` files
- 🚀 Fully automated *one-command* workflow

---

# 🚀 Quickstart

```bash
martinisurf \
    --pdb enzyme.pdb \
    --surface C1 --lx 20 --ly 20 \
    --ResA --ResB --dist 10

