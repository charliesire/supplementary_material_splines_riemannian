---
title: "Supplementary Material: Spline Interpolation on Riemannian Manifolds"
output: html_document
---

# supplementary_material_splines_riemannian

This repository contains the **Supplementary Material** for the article *"Spline Interpolation on Riemannian Manifolds"*.

## PDF Supplement

- **`supplementary_material.pdf`**  
  Presents the computation of the **mass and stiffness matrices** for 2-manifolds embedded in \(\mathbb{R}^3\), using natural local coordinate charts such as **cylindrical** or **spherical coordinates**.

## Code

The folder **`code_splines`** contains all the code associated with the paper *"Spline Interpolation on Compact Riemannian Manifolds"*. Its content includes:

- **`utils_harmonics.R`** – Functions for the method using **spherical harmonics**.  
- **`utils_matrix_triangu.R`** – Functions to build the different matrices associated with the triangulation, e.g., \(M\), \(F\), and \(A_n\).  
- **`utils_splines.R`** – Functions for **spline prediction** and **likelihood estimation** using the finite-element method.  
- **`utils_plot.R`** – Functions to **plot results in R**.  
- **`utils_plot_python.py`** – Functions to **plot results in Python**, especially for 3D visualizations.

### Sub-repositories

Two sub-folders implement the methods on specific manifolds:

1. **`cylinder`** – Codes for the method on the cylinder.  
2. **`sphere`** – Codes for the method on the sphere.  

Each of these sub-repositories contains:

- A **Jupyter notebook (`.ipynb`)** with the implementation.  
- **Results folders**: `first_scenario` and `second_scenario`.  
- **`plot_3d.ipynb`** – Notebook providing the 3D plots featured in the article.
