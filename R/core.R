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
  baseARI <- ARIs(merger$initialMat, unclustered = unclustered)
  # Normalize the impARI so that we take the mean over the same values.
  ARI <- merger$ImpARI *
    (ncol(merger$initialMat) - 1) / sum(upper.tri(baseARI))
  baseARI <- baseARI[upper.tri(baseARI)] %>% mean()
  ARI <- c(baseARI, ARI)
  ARI <- cumsum(ARI)
  return(ARI)
}

#' clusterConversion
#'
#' Find the conversion between the old cluster and the final clusters
#' @param merger the result from having run \code{\link{Dune}}
#'  on the dataset
#' @param p A value between 0 and 1. We stop when the mean ARI has improved by p
#' of the final total improvement. Default to 1 (i.e running the full merging).
#' @param n_steps Alternatively, you can specifiy the number of merging steps to
#' do before stopping.
#' @return A list containing a matrix per clustering method, with a column for the old labels
#' and a column for the new labels.
#' @importFrom magrittr %>%
#' @importFrom purrr walk
#' @importFrom dplyr n_distinct filter
#' @examples
#' data("clusMat", package = "Dune")
#' merger <- Dune(clusMat = clusMat)
#' clusterConversion(merger)[[2]]
#' @export
clusterConversion <- function(merger, p = 1, n_steps = NULL) {
  # Compute ARI imp and find where to stop the merge
  merges <- merger$merges
  if (is.null(n_steps)) {
    j <- whenToStop(merger, p = p)
  } else {
    j <- n_steps
  }
  merges <- merges[1:j, ]
  updates <- lapply(seq_len(ncol(merger$initialMat)), function(clusLab){
    clusters <- unique(merger$initialMat[, clusLab])
    update <- data.frame(old = clusters, new = clusters)
    return(update)
  })
  walk(seq_len(nrow(merges)), function(i) {
    clusLab <- merges[i, 1]
    clus <- max(merges[i, 2:3])
    id <- updates[[clusLab]]$new == clus
    updates[[clusLab]][id, "new"] <<- min(merges[i, 2:3])
  })
  return(updates)
}

#' Find the clustering matrix that we would get if we stopped the ARI merging
#' early
#' @param merger the result from having run \code{\link{Dune}}
#'  on the dataset
#' @param p A value between 0 and 1. We stop when the mean ARI has improved by p
#' of the final total improvement. Default to 1 (i.e running the full merging).
#' @param n_steps Alternatively, you can specifiy the number of merging steps to
#' do before stopping.
#' @return A matrix with the same dimensions as the currentmMat of the merger
#' argument
#' @examples
#' data("clusMat", package = "Dune")
#' merger <- Dune(clusMat = clusMat)
#' head(intermediateMat(merger, n_steps = 1))
#' @importFrom magrittr %>%
#' @import dplyr
#' @import tidyr
#' @export
intermediateMat <- function(merger, p = 1, n_steps = NULL) {
  # Get the conversion
  oldToNew <- clusterConversion(merger = merger, p = p, n_steps = n_steps)
  initialMat <- merger$initialMat %>% as.data.frame()
  names(oldToNew) <- colnames(initialMat)
  oldToNew <- bind_rows(oldToNew, .id = "cluster_label")
  initialMat <- initialMat %>%
    dplyr::mutate(cells = rownames(initialMat)) %>%
    pivot_longer(cols = colnames(initialMat), values_to = "old",
                 names_to = "cluster_label")
  newMat <- full_join(oldToNew, initialMat,
                      by = c("old" = "old", "cluster_label" = "cluster_label")) %>%
    select(cluster_label, new, cells) %>%
    pivot_wider(names_from = "cluster_label", values_from = "new") %>%
    arrange(cells)
  suppressWarnings(rownames(newMat) <- newMat$cells)
  newMat <- newMat %>% select(-cells)
  return(newMat)
}

#' Track the evolution of a function along merging
#'
#' For a given ARI merging, compute the evolution on the function f
#'
#' @param merger the result from having run \code{\link{Dune}}
#' on the dataset
#' @param f the function used. It must takes as input a clustering matrix and
#' return a value
#' @param p A value between 0 and 1. We stop when the mean ARI has improved by p
#' of the final total improvement. Default to 1 (i.e running the full merging).
#' @param n_steps Alternatively, you can specifiy the number of merging steps to
#' do before stopping.
#' @param ... additional arguments passed to f
#' @return a vector of length the number of merges
#' @importFrom magrittr %>%
#' @import dplyr
#' @examples
#' # Return the number of clusters for the fourth cluster label
#' data("clusMat", package = "Dune")
#' merger <- Dune(clusMat = clusMat)
#' f <- function(clusMat, i) dplyr::n_distinct(clusMat[, i])
#' functionTracking(merger, f, i = 4)
#' @export
functionTracking <- function(merger, f, p = 1, n_steps = NULL, ...){
  # Go over the merge and compute the function as we go
  if (is.null(n_steps)) {
    j <- whenToStop(merger, p = p)
  } else {
    j <- n_steps
  }
  values <- rep(0, j + 1)
  values[1] <- f(merger$initialMat, ...)
  for (i in seq_len(j)) {
    values[i + 1] <- f(intermediateMat(merger, n_steps = i), ...)
  }
  return(values)
}
