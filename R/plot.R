#' Plot an heatmap of the ARI matrix
#'
#' We can compute the ARI between pairs of cluster labels. This function plots
#' a matrix where a cell is the adjusted Rand Index between cluster label of
#' row i and cluster label of column j.
#' @param clusMat The clustering matrix with a row per cell and a column per
#' clustering label type
#' @param unclustered The value assigned to unclustered cells. Default to \code{NULL}
#' @param labels Whether to also display the ARI values. Default to TRUE
#' @return a \code{\link{ggplot}} object
#' @importFrom dplyr mutate
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @examples
#' data("clusMat", package = "Dune")
#' merger <- Dune(clusMat = clusMat)
#' plotARIs(merger$initalMat)
#' plotARIs(merger$currentMat)
#' @import ggplot2
#' @export

plotARIs <- function(clusMat, unclustered = NULL, labels = TRUE) {
  ARI <- ARIs(clusMat, unclustered = unclustered)
  p <- ARI %>% as.data.frame() %>%
    dplyr::mutate(label = rownames(ARI)) %>%
    tidyr::gather(key = label2, value = ari, -(ncol(ARI) + 1)) %>%
    ggplot(aes(x = label, y = label2, fill = ari)) +
    geom_tile() +
    scale_fill_viridis_c(limits = c(0, 1)) +
    theme_classic() +
    theme(axis.line = element_blank())
  if (labels) {
    p <- p  +
      geom_text(aes(label = round(ari, 2))) +
      guides(fill = FALSE)
  }
  return(p)
}

#' Plot the reduction in cluster size for an ARI merging with \code{Dune}
#' @param merger The output from an ARI merging, by calling \code{\link{Dune}}
#' @return a \code{\link{ggplot}} object
#' #' @importFrom dplyr mutate
#' @importFrom tidyr gather
#' @import ggplot2
#' @importFrom magrittr %>%
#' @examples
#' data("clusMat", package = "Dune")
#' merger <- Dune(clusMat = clusMat)
#' plotPrePost(merger)
#' @export
plotPrePost <- function(merger) {
  pre <- apply(merger$initalMat, 2, function(x) length(unique(x)))
  post <- apply(merger$currentMat, 2, function(x) length(unique(x)))
  df <- data.frame(methods = names(pre),
                   before = pre,
                   after = post) %>%
    tidyr::gather(key = "time", value = "Nb", -methods) %>%
    dplyr::mutate(time = factor(time, levels = c("before", "after")))
  p <- ggplot(df, aes(x = methods, y = Nb, fill = time)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_classic() +
    scale_fill_viridis_d(option = "E") +
    labs(x = "Clustering Methods", y = "Number of clusters", fill = "") +
    ggtitle("Reduction in number of clusters with ARI merging")
  return(p)
}

#' ARI improvement plot
#'
#' A plot to see how ARI improves over merging
#' @param merger the result from having run \code{\link{Dune}}
#'  on the dataset
#' @param unclustered The value assigned to unclustered cells. Default to
#' \code{NULL}
#' @return a \code{\link{ggplot}} object
#' @importFrom dplyr mutate
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @examples
#' data("clusMat", package = "Dune")
#' merger <- Dune(clusMat = clusMat)
#' ARItrend(merger)
#' @export
#' @import ggplot2
ARItrend <- function(merger, unclustered = NULL) {
  baseMat <- merger$initalMat
  ARI <- ARIImp(merger, unclustered = unclustered)
  n_clus <- lapply(1:nrow(merger$merges), function(m){
    diff <- rep(0, ncol(baseMat))
    diff[merger$merges[m, 1]] <- -1
    matrix(diff, nrow = 1)
  }) %>%
    do.call('rbind', args = .)
  n_clus <- rbind(sapply(baseMat, n_distinct) %>% matrix(data = ., nrow = 1),
                  n_clus)
  n_clus <- apply(n_clus, 2, cumsum)
  colnames(n_clus) <- colnames(baseMat)
  df <- data.frame(step = 0:length(merger$ImpARI),
                   ARI_Imp = ARI,
                   n_clus) %>%
    tidyr::gather(key = "change", value = "value", -step) %>%
    dplyr::mutate(type = ifelse(change == "ARI_Imp", "ARI Improvement",
                                "Nb of clusters"))
  p <- ggplot(df, aes(x = step, y = value)) +
    geom_path(size = 2, aes(group = change, col = change)) +
    facet_wrap(~type, scales = "free") +
    theme_classic() +
    scale_x_continuous(breaks = c(0, length(merger$ImpARI)),
                       labels = c("Initial", "Final")) +
    geom_hline(yintercept = min(ARI) + .9 * (max(ARI) - min(ARI)),
               col = "grey", linetype = "dashed", size = 2) +
    geom_vline(xintercept = min(which(ARI >= min(ARI) +
                                        .9 * (max(ARI) - min(ARI)))),
               col = "grey", linetype = "dashed", size = 2) +
    labs(y = "Change over merging",
         col = "type")
  return(p)
}
