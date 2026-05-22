Supplementary repository for Resin, Yang and Gneiting (2026+)
================
Johannes Resin

#### This repository provides supplementary code for the paper “Hierarchies of Calibration: Classification meets Regression” by Johannes Resin, Lu Yang and Tilmann Gneiting ([arXiv:XXX](https://arxiv.org/abs/XXX)).

##### Required packages

``` r
library(RcppAlgos)
library(pracma)
library(RColorBrewer)
library(tikzDevice)
library(matlib)
```

##### Load necessary functions

``` r
source("functions.R")
```

The repository implements the (theoretical) examples from the paper as
explained in Appendix B of the paper. The notebook ‘illustration.Rmd’
(and its rendered version ‘illustration.html’) provides an illustrative
example on how to use the functionality offered by the repository, where
we revisit Example 2.4 (b) from Gneiting and Resin (2023,
([DOI:10.1214/23-EJS2180](https://doi.org/10.1214/23-EJS2180))).
Properties of all other relevant examples from the paper are verified in
‘examples.Rmd’ (and the rendered version ‘examples.html’).
