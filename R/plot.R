#' Plot an heatmap of the ARI matrix
#'
#' We can compute the ARI between pairs of cluster labels. This function plots
#' a matrix where a cell is the adjusted Rand Index between cluster label of
#' row i and cluster label of column j.
#' @param ARI the matrix of pairwise ARI
#' @param small whether to also display the values
#' @return a \code{\link{ggplot}} object
#' @importFrom dplyr mutate
#' @importFrom tidyr gather
#' @import ggplot2
#' @export

plotARIs <- function(ARI, small = T) {
  p <- ARI %>% as.data.frame() %>%
    dplyr::mutate(label = rownames(ARI)) %>%
    tidyr::gather(key = label2, value = ari, -(ncol(ARI) + 1)) %>%
    ggplot(aes(x = label, y = label2, fill = ari)) +
    geom_tile() +
    scale_fill_viridis_c(limits = c(0, 1)) +
    theme_classic() +
    theme(axis.line = element_blank())
  if (small) {
    p <- p  +
      geom_text(aes(label = round(ari, 2))) +
      guides(fill = F)
  }
  return(p)
}

#' Plot the reduction in cluster size for an ARI merging
#' @param merger The output from an ARI merging
#' @return a \code{\link{ggplot}} object
#' #' @importFrom dplyr mutate
#' @importFrom tidyr gather
#' @import ggplot2
#' @export
plotPrePost <- function(merger) {
  r1 <- which(colnames(merger$initalMat) == "RsecT")
  if (all.equal(integer(0) ,r1) != TRUE) {
    pre <- apply(merger$initalMat[,-r1], 2, function(x) length(unique(x)))
  } else {
    pre <- apply(merger$initalMat, 2, function(x) length(unique(x)))
  }
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

#' plot the ARI improvement between methods
#'
#' The output from this function is a grid of 4 plots, based on the
#' \code{\link{plotARIs}} function. The first row of 2 plots is before the
#' merging procedure, the second is after the merging procedure. The first
#' column is when RSEC uses all assigned cells, the second column is when only
#' using the cells that RSEC cluster for computing the ARI between RSEC and
#' another partition of the data.
#' @param merger the result from having run \code{\link{mergeManyPairwise}}
#' on the dataset
#' @param Rsec default to TRUE. Whether there is an RSEC cluster labels(i.e one that does not cluster all cells)
#' @return the output from \code{\link{plot_grid}}
#' @importFrom dplyr mutate
#' @importFrom tidyr gather
#' @import ggplot2
#' @export
#' @importFrom cowplot plot_grid ggdraw draw_plot
plotARIReduce <- function(merger, RSEC = FALSE) {
  # Before, No unclustered cells for RSEC
  r1 <- which(colnames(merger$initalMat) == "RsecT")
  if (all.equal(integer(0) ,r1) != TRUE) {
    InitialARI <- apply(merger$initalMat[, -r1], 2, function(x) {
      apply(merger$initalMat[, -r1], 2, function(y) {
        inds <- x != -1 & y != -1
        xa <- x[inds]
        ya <- y[inds]
        adjustedRandIndex(xa, ya)
      })
    })

    p1 <- plotARIs(InitialARI) +
      ggtitle("ARI before any merging, no unclustered cells for RSEC") +
      theme(title = element_text(size = 6))

    # Before, All cells assigned
    r2 <- which(colnames(merger$initalMat) == "Rsec")
    InitialARI <- apply(merger$initalMat[,-r2], 2, function(x) {
      apply(merger$initalMat[,-r2], 2, function(y) {
        adjustedRandIndex(x, y)
      })
    })
    colnames(InitialARI)[colnames(InitialARI) == "RsecT"] <- "Rsec"
    rownames(InitialARI)[rownames(InitialARI) == "RsecT"] <- "Rsec"

    p2 <- plotARIs(InitialARI) +
      ggtitle("ARI before any merging, all cells assigned") +
      theme(title = element_text(size = 6))

    ## After, No unclustered cells for Rsec
    FinalARI <- apply(merger$currentMat, 2, function(x) {
      apply(merger$currentMat, 2, function(y) {
        inds <- x != -1 & y != -1
        xa <- x[inds]
        ya <- y[inds]
        adjustedRandIndex(xa, ya)
      })
    })

    p3 <- plotARIs(FinalARI) +
      ggtitle("ARI after merging, no unclustered cells for RSEC") +
      theme(title = element_text(size = 6))

    ## After, Rsec with all cells
    currentMat <- merger$currentMat
    currentMat[, "Rsec"] <- assignRsec(merger)
    FinalARI <- apply(currentMat, 2, function(x) {
      apply(currentMat, 2, function(y) {
        adjustedRandIndex(x, y)
      })
    })

    p4 <- plotARIs(FinalARI) +
      ggtitle("ARI after merging, all cells assigned") +
      theme(title = element_text(size = 6))

    p <- plot_grid(ggdraw() + draw_plot(p1),
              ggdraw() + draw_plot(p2),
              ggdraw() + draw_plot(p3),
              ggdraw() + draw_plot(p4),
              ncol = 2, rel_heights = rep(.25, 4), rel_widths = rep(.25, 4))
  } else {
    InitialARI <- apply(merger$initalMat, 2, function(x) {
      apply(merger$initalMat, 2, function(y) {
        adjustedRandIndex(x, y)
      })
    })

    p1 <- plotARIs(InitialARI) +
      ggtitle("ARI before any merging") +
      theme(title = element_text(size = 6))

    ## After
    FinalARI <- apply(merger$currentMat, 2, function(x) {
      apply(merger$currentMat, 2, function(y) {
        adjustedRandIndex(x, y)
      })
    })

    p2 <- plotARIs(FinalARI) +
      ggtitle("ARI after merging") +
      theme(title = element_text(size = 6))

    p <- plot_grid(ggdraw() + draw_plot(p1),
              ggdraw() + draw_plot(p2),
              ncol = 2, rel_heights = rep(1, 2), rel_widths = rep(.33, 2))
  }
  return(p)
}

#' ARI improvement plot
#'
#' A plot to see how ARI improves over merging
#' @param merger the result from having run \code{\link{mergeManyPairwise}}
#'  on the dataset
#' @return a \code{\link{ggplot}} object
#' @importFrom dplyr mutate
#' @importFrom tidyr gather
#' @export
#' @import ggplot2
ARItrend <- function(merger) {
  baseMat <- merger$initalMat
  j <- which(colnames(baseMat) == "RsecT")
  if (all.equal(integer(0) ,j) != TRUE) {
    baseMat <- baseMat[, -j]
  }
  ARI <- ARIImp(merger)
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

