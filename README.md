# Dune

<!-- badges: start -->
[![R-CMD-check](https://github.com/HectorRDB/Dune/workflows/R-CMD-check/badge.svg?branch=master)](https://github.com/HectorRDB/Dune/actions)
[![codecov](https://codecov.io/gh/HectorRDB/Dune/branch/master/graph/badge.svg?token=snxfXtj87B)](https://codecov.io/gh/HectorRDB/Dune)
<!-- badges: end -->

<p align="center">
  <img src="vignettes/logo.png" width="30%"/>
</p>

## Contents

- [Overview](#overview)
- [Installation](#installation)
- [Demo](#demo)
- [Issues](#Issues-and-bug-reports)
- [License](./LICENSE.md)
- [Citation](./inst/CITATION)


## Overview

Dune is an R Package that provides a parameter-free method for optimizing the trade-off between the resolutionof the clusters and their replicability across datasets. Dune  method takes as input a set of clustering results on a dataset, and iteratively merges clusters within those clusterings in order to maximize their concordance.  


## Installation

To install the current version of *Dune* , you will need Bioconductor 3.11 (i.e devel version):

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("Dune")
```
A version compatible with R 3.6 is also available with the following commands:

```
if(!requireNamespace("devtools", quietly = TRUE)) {
 install.packages("devtools") 
}
devtools::install_github("HectorRDB/Dune", ref = "before_bioc_release")
```

The installation should only take a few seconds.
The dependencies of the package are listed in the DESCRIPTION file of the package. 
To install R and Bioconductor, please refer to the [Bioconductor install page](https://www.bioconductor.org/install/).

## Issues and bug reports

Please use the [github issues](https://github.com/HectorRDB/Dune/issues) to submit issues, bug reports, and comments.

## Demo 

Start with the vignette [online](https://hectorRDB.github.io/Dune/articles/Dune.html).


If you to look at the vignette source code, you can either run `browseVignettes(package = "Dune")` if you install the package from Bioconductor find the code for the vignette [here](./vignettes/Dune.Rmd)
