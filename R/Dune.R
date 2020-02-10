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
Dune <- function(clusMat, unclustered = NULL, nCores = 1, verbose = FALSE) {
  # Initialize the values
  clusters <- lapply(as.data.frame(clusMat, stringsAsFactors = FALSE), unique)
  currentMat <- clusMat
  pairPartitions <- combn(colnames(clusMat), 2)
  confMats <- lapply(as.data.frame(pairPartitions, stringsAsFactors = FALSE),
                     function(pair){
    C1 <- pair[1]
    C2 <- pair[2]
    confusionMatrix <- table(clusMat[, C1], clusMat[, C2])
    ARI <- .adjustedRandIndex(confusionMatrix)
    return(list("C1" = C1, "C2" = C2, "confusionMatrix" = confusionMatrix,
                "ARI" = ARI))
  })

  working <- TRUE
  merges <- NULL
  ImpARI <- NULL

  # Try to see if any merge would increse
  while (working) {
    # For every cluster label list
    mergeResults <- parallel::mclapply(colnames(clusMat), function(C) {
      # Get all pairs of clusters
      clusterNames <- clusters[[C]]
      if (!is.null(unclustered)) {
        clusPairs <- combn(clusterNames[clusterNames != unclustered], 2)
      } else {
        clusPairs <- combn(clusterNames, 2)
      }

      # For every pair of labels in that list
      deltaARI <- apply(clusPairs, 2, function(pair) {
        m1 <- as.character(min(pair))
        m2 <- as.character(max(pair))
        localARIs <- lapply(confMats, function(confMat){
          if (confMat$C1 == C) {
            confusionMatrix <- confMat$confusionMatrix
            confusionMatrix[m1, ] <- confusionMatrix[m1, ] +
              confusionMatrix[m2, ]
            confusionMatrix <- confusionMatrix[
              -which(rownames(confusionMatrix) == m2), ]
            return(.adjustedRandIndex(confusionMatrix) - confMat$ARI)
          } else{
            if (confMat$C2 == C)  {
              confusionMatrix <- confMat$confusionMatrix
              confusionMatrix[, m1] <- confusionMatrix[, m1] +
                confusionMatrix[, m2]
              confusionMatrix <- confusionMatrix[,
                                      -which(colnames(confusionMatrix) == m2)]
              return(.adjustedRandIndex(confusionMatrix) - confMat$ARI)
            } else {
              return(0)
            }
          }
        })
        localARIs <- unlist(localARIs)
      })

      if (is.null(dim(deltaARI))) {
        return(deltaARI)
      } else {
        return(colMeans(deltaARI))
      }
    }, mc.cores = nCores)

    # Find best pair to merge
    maxs <- sapply(mergeResults, max)

    # Only merge if it improves ARI
    if (max(maxs) > 0) {
      # Find the cluster label where we merge
      clusLabel <- colnames(clusMat)[which.max(maxs)]
      # Find the two clusters to merge
      clusterNames <- clusters[[clusLabel]]
      if (!is.null(unclustered)) {
        clusPairs <- combn(clusterNames[clusterNames != unclustered], 2)
      } else {
        clusPairs <- combn(clusterNames, 2)
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
      m1 <- as.character(min(pair))
      m2 <- as.character(max(pair))
      confMats <- lapply(confMats, function(confMat){
        updatedConfMat <- confMat
        if (confMat$C1 == clusLabel) {
          confusionMatrix <- confMat$confusionMatrix
          confusionMatrix[m1, ] <- confusionMatrix[m1, ] +
            confusionMatrix[m2, ]
          confusionMatrix <- confusionMatrix[
            -which(rownames(confusionMatrix) == m2), ]
          updatedConfMat$confusionMatrix <- confusionMatrix
          updatedConfMat$ARI <- .adjustedRandIndex(confusionMatrix)
        } else{
          if (confMat$C2 == clusLabel)  {
            confusionMatrix <- confMat$confusionMatrix
            confusionMatrix[, m1] <- confusionMatrix[, m1] +
              confusionMatrix[, m2]
            confusionMatrix <- confusionMatrix[,
                                      -which(colnames(confusionMatrix) == m2)]
            updatedConfMat$confusionMatrix <- confusionMatrix
            updatedConfMat$ARI <- .adjustedRandIndex(confusionMatrix)
          }}
        return(updatedConfMat)
      })

      # tracking
      merges <- rbind(merges, c(clusLabel, pair))
      ImpARI <- c(ImpARI, max(maxs))
      if (verbose) {
        print(c(clusLabel, pair))
      }
    } else {
      working <- FALSE
    }

    # If no more to merge in any of them, stop
    if (sum(sapply(clusters, length) == 1) == length(clusters)) break()
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
  return(merger)
}
