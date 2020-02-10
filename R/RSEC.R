#' Assign cells using the assignUnassigned function of RSEC
#'
#' By default, RSEC does not assign all cells, leaving those as "-1".
#' However, to give a fair comparison with other labels, it is necessary to
#' assign cells so it is necessary to track that along the merging
#' @param merger the result from having run \code{\link{Dune}}
#' on the dataset
#' @param p when to stop the merging, when mean ARI has improved to p (between 0
#' and 1) of the final value.
#' @importFrom magrittr %>%
#' @return rsec are now all assigned
assignRsec <- function(merger, p = 1) {
  if (p < 0 | p > 1) {
    stop("p must be between zero and one")
  }
  if (p == 0) {
    return(merger$initialMat[,"RsecT"])
  }
  ARI <- ARIImp(merger)
  K <- min(which(ARI >= min(ARI) + p * (max(ARI) - min(ARI))))

  r1 <- which(colnames(merger$initialMat) == "RsecT")
  r2 <- which(colnames(merger$initialMat) == "Rsec")

  currentMat <- merger$currentMat
  Rsec_merges <- merger$merges[1:(K - 1), ]
  Rsec_merges <- Rsec_merges[Rsec_merges[,1] == 2, ]
  if (is.null(dim(Rsec_merges))) Rsec_merges <- matrix(Rsec_merges, nrow = 1)
  Rsec_merges <- Rsec_merges[, -1]
  if (is.null(dim(Rsec_merges))) Rsec_merges <- matrix(Rsec_merges, nrow = 1)
  if (nrow(Rsec_merges) == 0) {
    return(merger$initialMat$RsecT)
  } else {
    assign <- lapply(1:nrow(currentMat), function(i) {
      cell <- merger$initialMat[i, r2]
      cellT <- merger$initialMat[i, r1]
      for (j in 1:nrow(Rsec_merges)) {
        if (cellT %in% Rsec_merges[j, ]) cellT <- min(Rsec_merges[j, ])
      }
      return(cellT)
    }) %>%
      unlist()
    return(assign)
  }
}

#' Find the consensus clustering between three methods and return the consensus
#' @param clusMat The clustering matrix with a row per cell and a column per
#' clustering label type
#' @param large If the dataset is too large to be handled by \code{\link{makeConsensus}},
#'  we use a different method to get the consensus clustering. Default to FALSE.
#' @param ... Other arguments passed to \code{\link{makeConsensus}}
#' @return a vector of cluster assignations for the consensus.
#' @importFrom clusterExperiment makeConsensus
#' @importFrom magrittr %>%
#' @import dplyr
Consensus <- function(clusMat, large = FALSE, ...) {
  if (!large) {
    cellsConsensus <- suppressWarnings(
      clusterExperiment::makeConsensus(x = as.matrix(clusMat),
                                       clusterLabel = "makeConsensus",
                                       proportion = 2/3, minSize = 100, ...)
    )
    return(cellsConsensus$clustering)
  } else {
    collapse <- .dupRemove(t(clusMat))
    df <- clusterExperiment::makeConsensus(t(collapse$smallMat),
                                           proportion = 2/3, minSize = 0, ...)
    df <- df$clustering
    df <- data.frame(clusters = df, comb = 1:length(df))
    mapping <- data.frame(samples = 1:nrow(clusMat), comb = collapse$mapping)
    df <- df %>% inner_join(mapping) %>%
      arrange(samples) %>%
      group_by(clusters) %>%
      mutate(size = n()) %>%
      ungroup() %>%
      mutate(clusters = clusters * (size >= 500) + (-1) * (size < 500))
    return(df$clusters)
  }
}
