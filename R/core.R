#' Perform the ARI merging
#'
#' Compute the ARI between every pair of clustering labels after merging every possible pair of clusters. Find the one that improves the ARI merging the most, merge the pair. Repeat until there is no improvement.
#' @param clusMat the matrix of samples by clustering labels.
#' @param unclustered The value assigned to unclustered cells. Default to \code{NULL}
#' @param nCores number of cores to use when parallelizing. Relies on
#' \code{\link{mclapply}}. Default to 1.
#' @param verbose Whether or not the print cluster merging as it happens.
#' @return A list with four components: the initial matrix of clustering labels,
#'  the final matrix of clustering labels, the merge info matrix and the ARI
#'  improvement vector.
#' @details The Dune algorithm merges pairs of clusters in order to improve the
#' mean adjusted Rand Index with other clustering labels.
#' @seealso clusterConversion ARIImp
#' @importFrom parallel mclapply
#' @importFrom mclust adjustedRandIndex
#' @examples
#' data("clusMat", package = "Dune")
#' merger <- Dune(clusMat = clusMat)
#' # clusters 11 to 14 from cluster label 5 and 3 are subset of cluster 2 from
#' # other cluster labels. Designing cluster 2 as unclustered stops the merge of those.
#' merger2 <- Dune(clusMat = clusMat, unclustered = 2)
#' merger$merges
#' merger2$merges
#' @export
Dune <- function(clusMat, unclustered = NULL, nCores = 1,
                 verbose = FALSE) {
  # Turn the matrix into a numeric matrix
  # clusMat <- apply(clusteringMatrix, 2, function(x) {
  #   x[x != "-1"] <- as.numeric(factor(x[x != "-1"]))
  #   x[x == "-1"] <- -1
  #   x <- as.integer(x)
  # })

  # Initialize the values
  clusters <- lapply(seq_len(ncol(clusMat)), function(clus){
    unique(clusMat[, clus])
  })
  currentMat <- clusMat
  bestARI <- ARIs(currentMat, unclustered = unclustered)
  working <- TRUE
  merges <- NULL
  ImpARI <- NULL

  # Try to see if any merge would increse
  while (working) {
    # For every cluster label list
    mergeResults <- parallel::mclapply(seq_len(ncol(currentMat)),
                                       function(clusLabel) {
      clus <- currentMat[, clusLabel]
      clusterNames <- clusters[[clusLabel]]
      if (!is.null(unclustered)) {
        clusPairs <- combn(clusterNames[clusterNames != unclustered], 2)
      } else {
        clusPairs <- combn(clusterNames, 2)
      }

      # For every pair of labels in that list
      deltaARI <- apply(clusPairs, 2, function(pair) {
        sapply((1:ncol(currentMat))[-clusLabel], function(otherClus) {
          clus[clus %in% pair] <- max(clus) + 1
          # Compute the ARI once we merge the two clusters
          mclust::adjustedRandIndex(clus, currentMat[, otherClus])
        })
      }) - bestARI[clusLabel, -clusLabel]

      return(colMeans(deltaARI))
    }, mc.cores = nCores)

    # Find best pair to merge
    maxs <- sapply(mergeResults, max)

    # Only merge if it improves ARI
    if (max(maxs) > 0) {
      # Find the cluster label where we merge
      clusLabel <- which.max(maxs)
      # Find the two clusters to merge
      clusterNames <- clusters[[clusLabel]]
      if (!is.null(unclustered)) {
        clusPairs <- combn(clusterNames[clusterNames != unclustered], 2)
      } else {
        clusPairs <- combn(clusterNames, 2)
      }
      pair <- clusPairs[, which.max(mergeResults[[clusLabel]])]
      # Find which cells belonged to that cluster
      indsPair <- which(currentMat[, clusLabel] %in% pair)
      # Replace their cluster label with a new label
      currentMat[indsPair, clusLabel] <- min(pair)
      # Remove the old cluster from the list
      clusters[[clusLabel]] <-
        clusters[[clusLabel]][clusters[[clusLabel]] != max(pair)]

      # update bestARI
      bestARI <- ARIs(currentMat, unclustered = unclustered)

      # tracking
      merges <- rbind(merges, c(clusLabel, pair))
      ImpARI <- c(ImpARI, max(maxs))
      if (verbose) print(c(clusLabel, pair))
    } else {
      working <- FALSE
    }

    # If no more to merge in any of them, stop
    if (sum(sapply(clusters, length) == 1) == length(clusters)) break()
  }

  colnames(merges) <- c("clusteringLabel", "cluster1", "cluster2")
  return(list("initalMat" = clusMat,
              "currentMat" = as.data.frame(currentMat),
              "merges" = as.data.frame(merges),
              "ImpARI" = ImpARI))
}

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
#' @export
assignRsec <- function(merger, p = 1) {
  if (p < 0 | p > 1) {
    stop("p must be between zero and one")
  }
  if (p == 0) {
    return(merger$initalMat[,"RsecT"])
  }
  ARI <- ARIImp(merger)
  K <- min(which(ARI >= min(ARI) + p * (max(ARI) - min(ARI))))

  r1 <- which(colnames(merger$initalMat) == "RsecT")
  r2 <- which(colnames(merger$initalMat) == "Rsec")

  currentMat <- merger$currentMat
  Rsec_merges <- merger$merges[1:(K - 1), ]
  Rsec_merges <- Rsec_merges[Rsec_merges[,1] == 2, ]
  if (is.null(dim(Rsec_merges))) Rsec_merges <- matrix(Rsec_merges, nrow = 1)
  Rsec_merges <- Rsec_merges[, -1]
  if (is.null(dim(Rsec_merges))) Rsec_merges <- matrix(Rsec_merges, nrow = 1)
  if (nrow(Rsec_merges) == 0) {
    return(merger$initalMat$RsecT)
  } else {
    assign <- lapply(1:nrow(currentMat), function(i) {
      cell <- merger$initalMat[i, r2]
      cellT <- merger$initalMat[i, r1]
      for (j in 1:nrow(Rsec_merges)) {
        if (cellT %in% Rsec_merges[j, ]) cellT <- min(Rsec_merges[j, ])
      }
      return(cellT)
    }) %>%
      unlist()
    return(assign)
  }
}

#' ARI improvement
#'
#' Compute the ARI improvement over the ARI merging procedure
#' @param merger the result from having run \code{\link{Dune}}
#'  on the dataset
#' @param unclustered The value assigned to unclustered cells. Default to \code{NULL}
#' @return a vector with the mean ARI between methods at each merge
#' @seealso ARItrend
#' @importFrom magrittr %>%
#' @examples
#' data("clusMat", package = "Dune")
#' merger <- Dune(clusMat = clusMat)
#' plot(0:nrow(merger$merges), ARIImp(merger))
#' @export
ARIImp <- function(merger, unclustered = NULL) {
  baseMat <- merger$initalMat
  baseARI <- ARIs(baseMat, unclustered = unclustered)
  baseARI <- baseARI[upper.tri(baseARI)] %>% mean()
  ARI <- c(baseARI, merger$ImpARI)
  ARI <- cumsum(ARI)
  return(ARI)
}

#' clusterConversion
#'
#' Find the conversion between the old cluster and the final clusters
#' @param merger the result from having run \code{\link{Dune}}
#'  on the dataset
#' @param unclustered The value assigned to unclustered cells. Default to \code{NULL}
#' @return a vector with the mean ARI between methods at each merge
#' @importFrom magrittr %>%
#' @export
clusterConversion <- function(merger, unclustered = NULL) {
  baseMat <- merger$initalMat
  baseARI <- ARIs(baseMat, unclustered = unclustered)
  baseARI <- baseARI[upper.tri(baseARI)] %>% mean()
  ARI <- c(baseARI, merger$ImpARI)
  ARI <- cumsum(ARI)
  return(ARI)
}

#' Find the clustering matrix that we would get if we stopped the ARI merging
#' early
#' @param merger the result from having run \code{\link{Dune}}
#' on the dataset
#' @param p A value between 0 and 1. We stop when the mean ARI has improved by p
#' of the final total improvement
#' @return A matrix with the same dimensions as the currentmMat of the merger
#' argument
#' @importFrom magrittr %>%
#' @export
intermediateMat <- function(merger, p = .9) {
  # Compute ARI imp and find where to stop the merge
  ARI <- ARIImp(merger)
  int_merges <- merger$merges
  j <- min(which(ARI[2:length(ARI)] >= min(ARI) + p * (max(ARI) - min(ARI))))
  int_merges <- int_merges[1:j, ]
  assign <- sapply(colnames(merger$currentMat), function(clus) {
    J <- which(colnames(merger$initalMat) == clus)
    if (sum(int_merges[, 1] == J) == 0) {
      return(merger$initalMat[, clus])
    } else {
      clus_merges <- int_merges[int_merges[, 1] == J, ] %>%
        as.matrix() %>% matrix(ncol = 3)
      cells <- merger$initalMat[, clus]
      for (i in 1:nrow(clus_merges)) {
        indsPair <- which(cells %in% clus_merges[i, 2:3])
        cells[indsPair] <- min(clus_merges[i, 2:3])
      }
      return(cells)
    }
  })

  j <- which(colnames(merger$initalMat) == "RsecT")
  if (all.equal(integer(0) ,j) != TRUE) {
    colnames(assign) <- colnames(merger$initalMat)[-j]
  } else {
    colnames(assign) <- colnames(merger$initalMat)
  }

  return(assign)
}

#' Track the evolution of a function along merging
#'
#' For a given ARI merging, compute the evolution on the function f with
#' another partition
#' @param merger the result from having run \code{\link{Dune}}
#' on the dataset
#' @param f the function used, can be computed on any partition of the space
#' @param ... additional arguments passed to f
#' @return a matrix with a column per initial clustering, and a row per merge
#' with the f value computed
#' @importFrom magrittr %>%
#' @export
FTracking <- function(merger, f, ...){
  # Go over the merge and compute the homogeneity as we go
  baseMat <- merger$initalMat
  j <- which(colnames(baseMat) == "Rsec")
  if (all.equal(integer(0) ,j) != TRUE) {
    baseMat <- baseMat[, -j]
  }
  currentMat <- baseMat

  Evolution <- apply(baseMat, 2, f, ...) %>%  matrix(ncol = ncol(baseMat))

  for (m in seq_len(nrow(merger$merges))) {
    wClus <- merger$merges[m, 1]
    clus <- currentMat[, wClus]
    pair <- merger$merges[m, 2:3]
    clus[clus %in% pair] <- min(pair)
    currentMat[, wClus] <- clus
    Evolution <- rbind(Evolution, Evolution[nrow(Evolution), ])
    Evolution[nrow(Evolution), wClus] <- f(clus, ...)
  }
  return(Evolution)
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
#' @export
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
