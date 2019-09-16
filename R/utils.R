#' Remove the duplicated cells which belong to the same cluster over all labels.
#' @param clusMat The clustering matrix with a row per cell and a column per
#' clustering label type
#' @return a list containing the reduced matrix and a way to map back to the original
.dupRemove <- function(clusMat) {
  whDup <- which(duplicated(t(clusMat)))
  val <- apply(clusMat, 2, paste, collapse = ",") # all combinations
  if (length(whDup) > 0) {
    clusMat <- clusMat[, -whDup, drop = FALSE]
  }
  valSm <- apply(clusMat, 2, paste, collapse = ",") # unique combinations
  ind <- match(val, valSm)
  return(list(smallMat = clusMat, mapping = ind))
}

#' @title ARI Matrix
#' @param clusMat The clustering matrix with a row per cell and a column per
#' clustering label type
#' @param unclustered The value assigned to unclustered cells. Default to \code{NULL}
#' @return The ARI matrix
#' @details In the ARI matrix where each cell **i,j** is the adjusted Rand Index
#' between columns **i** and **j** of the original \code{clusMat}.
#' If \code{unclustered} is not NULL, the cells which have been assigned to the
#' \code{unclustered} cluster will not be counted towards computing the ARI.
#' @importFrom mclust adjustedRandIndex
#' @export
ARIs <- function(clusMat, unclustered = NULL) {
  ARI <- apply(clusMat, 2, function(x) {
    apply(clusMat, 2, function(y) {
      if (!is.null(unclustered)) {
        x_unc <- x[x != unclustered & y != unclustered]
        y_unc <- y[x != unclustered & y != unclustered]
      } else {
        x_unc <- x
        y_unc <- y
      }
      mclust::adjustedRandIndex(x_unc, y_unc)
    })
  })
  if (is.null(colnames(clusMat))) {
    rownames(ARI) <- colnames(ARI) <- 1:ncol(clusMat)
  } else {
    rownames(ARI) <- colnames(ARI) <- colnames(clusMat)
  }

  return(ARI)
}

#' @title ARI Matrix
#' @param merger the result from having run \code{\link{Dune}}
#'  on the dataset
#' @param p A value between 0 and 1. We stop when the mean ARI has improved by p
#' of the final total improvement. Default to 1 (i.e running the full merging).
#' @return An integer giving the step where to stop.
#' @details The \code{\link{Dune}} process improves the mean ARI. This return
#' the first merging step after which the mean ARI has been improved by p of the
#' total. Setting p = 1 just return the number of merges.
#' @importFrom mclust adjustedRandIndex
#' @examples
#' data("clusMat", package = "Dune")
#' merger <- Dune(clusMat = clusMat)
#' whenToStop(merger, p = .5)
#' @export
whenToStop <- function(merger, p) {
  if (p < 0 | p > 1) stop("p must be between 0 and 1")
  ARI <- ARIImp(merger)
  j <- min(which(ARI[2:length(ARI)] >= min(ARI) + p * (max(ARI) - min(ARI))))
  return(j)
}
