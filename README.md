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
- 🧱 **Integration with martinize2 and Gō–Martini**
- 🧰 **Full pipeline mode:**  
  Generates *all* topology, index, restraints, and MDP files


## 🔧 **Installation**

- conda create -n martinisurf
- conda activate martinisurf
- pip install -r surfmartini
- conda install pip
- pip install .

### 🔍 What this does

- 📥 Downloads **1RJW** from RCSB  
- 🧬 Runs **martinize2** to generate the coarse-grained model  
- 🧱 Builds a **20 × 20 nm** Martini surface  
- 📐 Orients the enzyme using **two anchor groups**  
- 🧰 Generates complete **GROMACS-ready simulation files**  

