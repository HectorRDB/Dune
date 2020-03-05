.findMergeResults <- function(C, clusters, unclustered, confMats) {
  # Get all pairs of clusters
  clusterNames <- clusters[[C]]
  if (!is.null(unclustered)) {
    if (length(clusterNames[clusterNames != unclustered]) == 1) {
      return(rep(0, length(confMats)))
    }
    clusPairs <- utils::combn(clusterNames[clusterNames != unclustered], 2)
  } else {
    if (length(clusterNames) == 1) {
      return(rep(0, length(confMats)))
    }
    clusPairs <- utils::combn(clusterNames, 2)
  }

  # For every pair of cluster in that list, compute how the ARI change if
  # we merge
  deltaARI <- apply(clusPairs, 2, .localARI, confMats = confMats, C = C)

  # Needed for the cases with only two partitions
  if (is.null(dim(deltaARI))) {
    return(deltaARI)
  } else {
    return(colMeans(deltaARI))
  }
}

.Dune <- function(clusMat,
                  unclustered = NULL,
                  verbose = FALSE,
                  parallel = FALSE,
                  BPPARAM = BiocParallel::bpparam()){
  # Initialize the values
  ## Unique cluster labels for all partitions
  currentMat <- clusMat
  if (is.null(colnames(clusMat))) colnames(currentMat) <- seq_len(ncol(clusMat))
  clusters <- lapply(as.data.frame(currentMat, stringsAsFactors = FALSE),
                     unique)
  pairPartitions <- utils::combn(colnames(currentMat), 2)
  # All confusion matrices with associated partitions and ARIs
  confMats <- .confusionMatrices(pairPartitions, currentMat)

  # This is used to keep track of all the info
  working <- TRUE
  merges <- NULL
  ImpARI <- NULL

  # Try to see if any merge would increse
  while (working) {
    # For every partition
    if (parallel) {
      mergeResults <- BiocParallel::bplapply(
        colnames(currentMat), .findMergeResults, clusters = clusters,
        unclustered = unclustered, confMats = confMats, BPPARAM = BPPARAM
      )
    } else {
      mergeResults <- lapply(colnames(currentMat), .findMergeResults,
                             clusters = clusters, unclustered = unclustered,
                             confMats = confMats)
    }


    # Find best pair to merge
    maxs <- sapply(mergeResults, max)

    # Only merge if it improves ARI
    if (max(maxs) > 0) {
      # Find the partition where we merge
      clusLabel <- colnames(currentMat)[which.max(maxs)]
      # Find the two clusters to merge
      clusterNames <- clusters[[clusLabel]]
      if (!is.null(unclustered)) {
        clusPairs <- utils::combn(clusterNames[clusterNames != unclustered], 2)
      } else {
        clusPairs <- utils::combn(clusterNames, 2)
      }
      pair <- clusPairs[, which.max(mergeResults[[which.max(maxs)]])]
      # Find which cells belonged to that cluster
      indsPair <- which(currentMat[, clusLabel] %in% pair)
      # Replace their cluster label with a new label
      currentMat[indsPair, clusLabel] <- min(pair)
      # Remove the old cluster from the list
      clusters[[clusLabel]] <-
        clusters[[clusLabel]][clusters[[clusLabel]] != max(pair)]

      # Update the confusion matrices
      confMats <- .updatedConfMats(pair, confMats, clusLabel)
      # Tracking
      merges <- rbind(merges, c(clusLabel, pair))
      ImpARI <- c(ImpARI, max(maxs))
      if (verbose) {
        print(c(clusLabel, pair))
      }
    } else {
      working <- FALSE
    }

    # If no more to merge in any of them, stop
    if (any(sapply(clusters, length) == 1)) {
      if (verbose) {
        print("We merge one of the partition entirely")
      }
      break()
    }
  }

  if (is.null(merges)) {
    stop("This resulted in no merges. Check input cluster labels")
  }
  colnames(merges) <- c("clusteringLabel", "cluster1", "cluster2")
  merger <- list(
    "initialMat" = clusMat,
    "currentMat" = as.data.frame(currentMat, stringsAsFactors = FALSE),
    "merges" = as.data.frame(merges, stringsAsFactors = FALSE),
    "ImpARI" = ImpARI
  )
  if (is.null(colnames(clusMat))) {
    colnames(merger$initialMat) <- colnames(currentMat)
  }
  return(merger)
}

#' Perform the ARI merging
#'
#' Compute the ARI between every pair of clustering labels after merging every possible pair of clusters. Find the one that improves the ARI merging the most, merge the pair. Repeat until there is no improvement.
#' @param clusMat the matrix of samples by clustering labels.
#' @param cluster_columns if \code{clusMat} is a \code{\link{SummarizedExperiment}},
#'  then this defines the columns of \code{colData} that are outputs from a clustering algorithm.
#' @param unclustered The value assigned to unclustered cells. Default to \code{NULL}
#' @param parallel Logical, defaults to FALSE.
#' Set to TRUE if you want to parallellize the fitting.
#' @param BPPARAM object of class \code{bpparamClass} that specifies the
#'   back-end to be used for computations.
#'   See \code{bpparam} in \code{BiocParallel} package for details.
#'   Won't be used if \code{parallel} is FALSE.
#' @param verbose Whether or not the print cluster merging as it happens.
#' @return A list with four components: the initial matrix of clustering labels,
#'  the final matrix of clustering labels, the merge info matrix and the ARI
#'  improvement vector.
#' @details The Dune algorithm merges pairs of clusters in order to improve the
#' mean adjusted Rand Index with other clustering labels. It returns a list with
#'  four components.:
#'  #' \itemize{
#'   \item \code{initialMat}: The initial matrix of cluster labels
#'   \item \code{currentMat}: The final matrix of cluster labels
#'   \item \code{merges}: The step-by-step detail of the merges, recapitulating
#'   which clusters where merged in which cluster label
#'   \item \code{impARI}: How much each merge improved the mean ARI between the
#'   cluster label that has been merged and the other cluster labels.
#' }
#' @seealso clusterConversion ARIImp
#' @importFrom BiocParallel bplapply bpparam
#' @import utils
#' @examples
#' data("clusMat", package = "Dune")
#' merger <- Dune(clusMat = clusMat)
#' # clusters 11 to 14 from cluster label 5 and 3 are subset of cluster 2 from
#' # other cluster labels. Designing cluster 2 as unclustered therefore means we
#' # do fewer merges.
#' merger2 <- Dune(clusMat = clusMat, unclustered = 2)
#' merger$merges
#' merger2$merges
#' @export
#' @rdname Dune
setMethod(f = "Dune",
          signature = c(clusMat = "matrix"),
          definition = function(clusMat,
                                unclustered = NULL,
                                verbose = FALSE,
                                parallel = FALSE,
                                BPPARAM = BiocParallel::bpparam()) {
            merger <- .Dune(clusMat = clusMat, unclustered = unclustered,
                            verbose = verbose, BPPARAM = BPPARAM)
            return(merger)
})

#' @rdname Dune
#' @importFrom BiocParallel bplapply bpparam
setMethod(f = "Dune",
          signature = c(clusMat = "data.frame"),
          definition = function(clusMat,
                                unclustered = NULL,
                                verbose = FALSE,
                                parallel = FALSE,
                                BPPARAM = BiocParallel::bpparam()) {
            merger <- .Dune(clusMat = clusMat, unclustered = unclustered,
                            verbose = verbose, BPPARAM = BPPARAM)
            return(merger)
          }
)

#' @rdname Dune
#' @importFrom SummarizedExperiment colData
#' @importFrom BiocParallel bplapply bpparam
setMethod(f = "Dune",
          signature = c(clusMat = "SummarizedExperiment"),
          definition = function(clusMat,
                                cluster_columns,
                                unclustered = NULL,
                                verbose = FALSE,
                                parallel = FALSE,
                                BPPARAM = BiocParallel::bpparam()) {
            df <- colData(clusMat)
            if (any(!cluster_columns %in% colnames(df))) {
              stop("All elements of cluster_columns should be in colData(clusMat)")
            }
            df <- as.data.frame(df)
            df <- df[,cluster_columns]
            merger <- Dune(clusMat = df, unclustered = unclustered,
                           verbose = verbose, BPPARAM = BPPARAM)
            return(merger)
          }
)
