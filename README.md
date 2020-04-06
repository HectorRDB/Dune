# Dune

<!-- badges: start -->
  [![Build Status](https://travis-ci.com/HectorRDB/Dune.svg?token=zVqScmB6wvrS9ZJSK57p&branch=before_bioc_release)](https://travis-ci.com/HectorRDB/Dune) [![codecov](https://codecov.io/gh/HectorRDB/Dune/branch/master/graph/badge.svg?token=snxfXtj87B)](https://codecov.io/gh/HectorRDB/Dune)
<!-- badges: end -->

<p align="center">
  <img src="vignettes/logo.png" width="30%"/>
</p>

## Dune overview

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

## Issues and bug reports

Please use the [github issues](https://github.com/HectorRDB/Dune/issues) to submit issues, bug reports, and comments.

## Usage 

Start with the vignette [online](https://hectorRDB.github.io/Dune/articles/Dune.html).
