#' Perform the ARI merging
#'
#' Compute the ARI between every pair of clustering labels after merging every possible pair of clusters. Find the one that improves the ARI merging the most, merge the pair. Repeat until there is no improvement.
#' @param clusteringMatrix the matrix of samples by clustering labels.
#' @param nCores number of cores to use when parallelizing. Relies on \code{\link{mclapply}}.
#' @return a list with the initial matrix of clustering labels, the final matrix of clustering labels, the merge info matrix and the ARI improvement.
#' @importFrom parallel mclapply
#' @importFrom mclust adjustedRandIndex
#' @export
mergeManyPairwise <- function(clusteringMatrix, nCores = 3) {
  # Turn the matrix into a numeric matrix
  clusMat <- apply(clusteringMatrix, 2, function(x) {
    x[x != "-1"] <- as.numeric(factor(x[x != "-1"]))
    x[x == "-1"] <- -1
    x <- as.integer(x)
  })

  # Initialize the values
  clusters <- apply(clusMat, 2, unique)
  rownames(clusMat) <- rownames(clusteringMatrix)
  currentMat <- clusMat
  baseARI <- apply(clusMat, 2, function(x) {
    apply(clusMat, 2, function(y) {
      mclust::adjustedRandIndex(x, y)
    })
  })
  bestARI <- baseARI
  working <- TRUE
  merges <- NULL
  ImpARI <- NULL

  # Try to see if any merge would increse
  while (working) {
    # Test all pairwise clusters to merge
    # For every cluster label list
    mergeResults <- parallel::mclapply(1:ncol(currentMat), function(whClus) {
      clus <- currentMat[, whClus]
      clusternames <- clusters[[whClus]]
      clusPairs <- combn(clusternames[clusternames != -1], 2)

      # For every pair of labels in that list
      deltaARI <- apply(clusPairs, 2, function(pair) {
        sapply((1:ncol(clusMat))[-whClus], function(otherClus) {
          clus[clus %in% pair] <- max(clus) + 1
          mclust::adjustedRandIndex(clus, currentMat[, otherClus])
        })
      }) - bestARI[whClus, -whClus]

      return(colMeans(deltaARI))
    }, mc.cores = nCores)

    # Find best pair to merge
    maxs <- sapply(mergeResults, max)

    # Only merge if it improves ARI
    if (max(maxs) > 0) {
      whClus <- which.max(maxs)
      # update clusters
      clusternames <- clusters[[whClus]]
      clusPairs <- combn(clusternames[clusternames != -1], 2)
      pair <- clusPairs[, which.max(mergeResults[[whClus]])]
      indsPair <- which(currentMat[, whClus] %in% pair)
      currentMat[indsPair, whClus] <- min(pair)
      clusters[[whClus]] <- unique(currentMat[, whClus])

      # update bestARI
      newARIs <- sapply((1:ncol(currentMat))[-whClus], function(j) {
        mclust::adjustedRandIndex(currentMat[, whClus], currentMat[, j])
      })
      bestARI[whClus, -whClus] <- newARIs
      bestARI[-whClus, whClus] <- newARIs

      # tracking
      merges <- rbind(merges, c(whClus, pair))
      ImpARI <- c(ImpARI, max(maxs))
      print(c(whClus, pair))
    } else {
      working <- FALSE
    }

    # If no more to merge in any of them, stop
    if (sum(sapply(clusters, length) == 1) == length(clusters)) stop()
  }

  colnames(merges) <- c("clustering", "cluster1", "cluster2")
  return(list("initalMat" = clusteringMatrix,
              "currentMat" = as.data.frame(currentMat),
              "merges" = merges,
              "ImpARI" = ImpARI))
}

#' Assign cells using the assignUnassigned function of RSEC
#'
#' By default, RSEC does not assign all cells, leaving those as "-1".
#' However, to give a fair comparison with other labels, it is necessary to
#' assign cells so it is necessary to track that along the merging
#' @param merger the result from having run \code{\link{mergeManyPairwise}}
#' on the dataset
#' @param p when to stop the merging, when mean ARI has improved to p (between 0
#' and 1) of the final value.
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
#' @param merger the result from having run \code{\link{mergeManyPairwise}}
#'  on the dataset
#' @return a vector with the mean ARI between methods at each step
#' @export
ARIImp <- function(merger) {
  baseMat <- merger$initalMat
  j <- which(colnames(baseMat) == "RsecT")
  if (all.equal(integer(0) ,j) != TRUE) {
    baseMat <- baseMat[, -j]
  }
  baseARI <- apply(baseMat, 2, function(x) {
    apply(baseMat, 2, function(y) {
      adjustedRandIndex(x, y)
    })
  })
  baseARI <- baseARI[upper.tri(baseARI)] %>% mean()
  ARI <- c(baseARI, merger$ImpARI)
  ARI <- cumsum(ARI)
  return(ARI)
}

#' Find the clustering matrix that we would get if we stopped the ARI merging
#' early
#' @param merger the result from having run \code{\link{mergeManyPairwise}}
#' on the dataset
#' @param p A value between 0 and 1. We stop when the mean ARI has improved by p
#' of the final total improvement
#' @return A matrix with the same dimensions as the currentmMat of the merger
#' argument
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
#' @param merger the result from having run \code{\link{mergeManyPairwise}}
#' on the dataset
#' @param f the function used, can be computed on any partition of the space
#' @param ... additional arguments passed to f
#' @return a matrix with a column per initial clustering, and a row per merge
#' with the f value computed
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
#' @param large If the dataset is too large to be handled by \code{\link{MakeConsensus}},
#'  we use a different method to get the consensus clustering. Default to FALSE.
#' @param ... Other arguments passed to \code{\link{MakeConsensus}}
#' @return a vector of cluster assignations for the consensus.
#' @importFrom clusterExperiment makeConsensus
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

#' Remove the duplicated cells which belong to the same cluster over all labels.
#' @param clusMat The clustering matrix with a row per cell and a column per
#' clustering label type
#' @return a list containing the reduced matri and a way to map back to the original
.dupRemove <- function(clusterMat) {
  whDup <- which(duplicated(t(clusterMat)))
  val <- apply(clusterMat, 2, paste, collapse = ",") # all combinations
  if (length(whDup) > 0) {
    clusterMat <- clusterMat[, -whDup, drop = FALSE]
  }
  valSm <- apply(clusterMat, 2, paste, collapse = ",") # unique combinations
  ind <- match(val, valSm)
  return(list(smallMat = clusterMat, mapping = ind))
}
