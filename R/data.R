#' A clustering matrix used to demonstrate the ari-merging process.
#'
#' @details This matrix has 100 samples with 5 cluster labels. Cluster labels 2
#'  trought 5 are modified versions of cluster label 1, where some clusters from
#'  label 1 where broken down into smaller clusters. 
#'  It is just a toy dataset that can be re-generated with the code in
#'  https://github.com/HectorRDB/Pipeline_Brain/blob/master/Sandbox/createToyDataset.R
#'
"clusMat"

#' Cluster labels for a subset of the allen Smart-Seq nuclei dataset
#'
#' @details This matrix of clusters was obtained by running 3 clustering algorithms
#' on a brain snRNA-Seq dataset from Tasic et .al (https://doi.org/10.1038/s41586-018-0654-5).
#' This dataset was then subsetted to the GABAergic neurons. 
#' Code to reproduce all this can be found in the github repository from the 
#' Dune paper (https://github.com/HectorRDB/Dune_Paper).
#'
"nuclei"
