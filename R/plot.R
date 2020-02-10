#' Plot an heatmap of the ARI matrix
#'
#' We can compute the ARI between pairs of cluster labels. This function plots
#' a matrix where a cell is the adjusted Rand Index between cluster label of
#' row i and cluster label of column j.
#' @param clusMat The clustering matrix with a row per cell and a column per
#' clustering label type
#' @param unclustered The value assigned to unclustered cells. Default to \code{NULL}
#' @param values Whether to also display the ARI values. Default to TRUE.
#' @param numericalLabels Whether labels are numerical values. Default to FALSE.
#' @return a \code{\link{ggplot}} object
#' @importFrom dplyr mutate
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @examples
#' data("clusMat", package = "Dune")
#' merger <- Dune(clusMat = clusMat)
#' plotARIs(merger$initialMat)
#' plotARIs(merger$currentMat)
#' @import ggplot2
#' @import RColorBrewer
#' @export

plotARIs <- function(clusMat, unclustered = NULL, values = TRUE,
                     numericalLabels = FALSE) {
  ARI <- ARIs(clusMat, unclustered = unclustered)
  df <- ARI %>% as.data.frame() %>%
    dplyr::mutate(label1 = rownames(ARI)) %>%
    tidyr::gather(key = label2, value = ari, -(ncol(ARI) + 1))
  if (numericalLabels) {
    df <- df %>%
      dplyr::mutate(label1 = as.numeric(label1), label2 = as.numeric(label2))
  }
  p <- ggplot(df, aes(x = label1, y = label2, fill = ari)) +
    geom_tile() +
    scale_fill_gradientn(colours = brewer.pal(9, "Spectral"),
                         limits = c(0, 1)) +
    theme_classic() +
    theme(axis.line = element_blank())
  if (values) {
    p <- p  +
      geom_text(aes(label = round(ari, 2)), size = 4) +
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
  pre <- apply(merger$initialMat, 2, function(x) length(unique(x)))
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
  baseMat <- merger$initialMat
  ARI <- ARIImp(merger, unclustered = unclustered)
  n_clus <- lapply(seq_len(nrow(merger$merges)), function(m){
    diff <- rep(0, ncol(baseMat))
    names(diff) <- colnames(merger$initialMat)
    diff[merger$merges[m, 1]] <- -1
    matrix(diff, nrow = 1)
  }) %>%
    do.call('rbind', args = .)
  n_clus <- rbind(apply(baseMat, 2, n_distinct) %>% matrix(data = ., nrow = 1),
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
    scale_color_brewer(type = "qual") +
    scale_x_continuous(breaks = c(0, length(merger$ImpARI)),
                       labels = c("Initial", "Final")) +
    labs(y = "Change over merging",
         col = "Type")
  return(p)
}

#' Plot confusion matrix
#'
#' A plot to visualize how alike two clustering labels are
#' @param x A vector of clustering labels or a matrix of clustering labels. See details.
#' @param y Optional. Another vector of clustering labels
#' @return a \code{\link{ggplot}} object
#' @importFrom dplyr mutate group_by
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @importFrom RColorBrewer brewer.pal
#' @examples
#' data("nuclei", package = "Dune")
#' ConfusionPlot(nuclei[, c("SC3", "Monocle")])
#' @export
#' @import ggplot2

ConfusionPlot <- function(x, y = NULL) {
  if (is.null(y)) {
    y <- x[, 2]
    x <- x[, 1]
  } else {
    if (length(x) != length(y)) {
      stop("x and y must have the same length")
    }
  }
  df <- table(x, y) %>%
    as.data.frame() %>%
    group_by(x) %>%
    mutate(total_x = sum(Freq)) %>%
    group_by(y) %>%
    mutate(total_y = sum(Freq),
           union = total_x + total_y - Freq,
           overlap = Freq / union) %>%
    ungroup() %>%
    arrange(desc(Freq)) %>%
    filter(Freq > 0)
  df$x <- factor(df$x, levels = unique(df$x))
  df$y <- factor(df$y, levels = unique(df$y))
  p <- ggplot(df, aes(x = x, y = y, col = overlap, size = Freq)) +
    geom_point() +
    theme_bw() +
    theme(legend.position = "top",
          rect = element_blank(),
          panel.border = element_blank(),
          legend.box.spacing = unit(0, units = "npc"),
          legend.margin	=  margin(r = .1, l = .1, unit = "npc")) +
    NULL +
    labs(col = "% of Overlap", size = "# of Cells") +
    scale_color_gradientn(colours = brewer.pal(11, "Spectral")) +
    guides(size = guide_legend(title.position = "top", fill = "grey"),
           col = guide_colourbar(title.position = "top",
                                 barwidth = unit(.2, "npc")))
  return(p)
  }


#' Plot the evolution of the pairwise ARIs as merging
#' happens
#'
#' Animated version of \code{\link{plotARIs}}
#'
#' @param merger the result from having run \code{\link{Dune}}
#'  on the dataset
#' @param unclustered The value assigned to unclustered cells. Default to
#' \code{NULL}
#' @param values Whether to also display the ARI values. Default to TRUE.
#' @param numericalLabels Whether labels are numerical values. Default to FALSE.
#' @param state_length Time between steps. Default to 1. See \code{\link{transition_states}}
#' for details.
#' @return a \code{gganim} object
#' @details See \code{\link{plotARIs}} and \code{\link{animate}}.
#' @importFrom purrr map
#' @importFrom tidyr gather
#' @importFrom dplyr mutate n_distinct
#' @importFrom magrittr %>%
#' @importFrom RColorBrewer brewer.pal
#' @importFrom gganimate transition_states
#' @import ggplot2
#' @examples
#' \dontrun{
#'   data("clusMat", package = "Dune")
#'   merger <- Dune(clusMat = clusMat)
#'   ARIEvolution(merger)}
#' @export
ARIEvolution <- function(merger, unclustered = NULL, values = TRUE,
                            numericalLabels = FALSE, state_length = 1) {
  ARI_matrices <- purrr::map_df(0:length(merger$ImpARI), function(step){
    ARI <- ARIs(clusMat = intermediateMat(merger, n_steps = step),
         unclustered = unclustered)
    ARI <- ARI %>%
      as.data.frame() %>%
      dplyr::mutate(label1 = rownames(ARI)) %>%
      tidyr::gather(key = label2, value = ari, -(ncol(ARI) + 1)) %>%
      dplyr::mutate(step = step)
  })
  if (numericalLabels) {
    ARI_matrices <- ARI_matrices %>%
      dplyr::mutate(label1 = as.numeric(label1), label2 = as.numeric(label2))
  }
  ARI_matrices <- ARI_matrices %>%
    dplyr::arrange(step)
  p <- ggplot(ARI_matrices, aes(x = label1, y = label2, fill = ari)) +
    geom_tile() +
    scale_fill_gradientn(colours = RColorBrewer::brewer.pal(9, "Spectral"),
                         limits = c(0, 1)) +
    theme_classic() +
    theme(axis.line = element_blank())
  if (values) {
    p <- p  +
      geom_text(aes(label = round(ari, 2)), size = 4) +
      guides(fill = FALSE)
  }
  p <- p +
  gganimate::transition_states(step,
                               transition_length = 0,
                               state_length =  state_length /
                                 dplyr::n_distinct(ARI_matrices$label1)) +
    ggtitle(paste0('Step {closest_state} of ', max(ARI_matrices$step)))
  return(p)
}

#' Plot the evolution of the ConfusionPlot as merging
#' happens
#'
#' Animated version of \code{\link{ConfusionPlot}}
#' @param merger the result from having run \code{\link{Dune}}
#'  on the dataset
#' @param unclustered The value assigned to unclustered cells. Default to \code{NULL}
#' @param x The name of the first cluster label to plot
#' @param y The name of the second cluster label to plot
#' @param state_length Time between steps. Default to 1. See \code{\link{transition_states}}
#' for details.
#' @return a \code{gganim} object
#' @details See \code{\link{ConfusionPlot}} and \code{\link{animate}}.
#' @importFrom purrr map
#' @import tidyr dplyr
#' @importFrom magrittr %>%
#' @importFrom RColorBrewer brewer.pal
#' @import ggplot2
#' @importFrom gganimate transition_states
#' @examples
#' \dontrun{
#'   data("clusMat", package = "Dune")
#'   merger <- Dune(clusMat = clusMat)
#'   ConfusionEvolution(merger, x = "A", y = "B")}
#' @export
ConfusionEvolution <- function(merger, unclustered = NULL, x, y, state_length = 1) {
  Freqs <- purrr::map_df(0:length(merger$ImpARI), function(step){
    clusMat <- intermediateMat(merger, n_steps = step) %>%
      as.matrix()
    df <- table(x = clusMat[, x], y = clusMat[, y]) %>%
      as.data.frame() %>%
      mutate(x = as.numeric(x),
             y = as.numeric(y)) %>%
      group_by(x) %>%
      mutate(total_x = sum(Freq),
             step = step) %>%
      group_by(y) %>%
      mutate(total_y = sum(Freq),
             union = total_x + total_y - Freq,
             overlap = Freq / union) %>%
      ungroup() %>%
      arrange(desc(Freq)) %>%
      filter(Freq > 0)
  })
  Freqs$x <- factor(Freqs$x, levels = sort(unique(Freqs$x)))
  Freqs$y <- factor(Freqs$y, levels = sort(unique(Freqs$y)))

  p <- ggplot(Freqs, aes(x = x, y = y, col = overlap, size = Freq)) +
    geom_point() +
    theme_bw() +
    theme(legend.position = "top",
          rect = element_blank(),
          panel.border = element_blank(),
          legend.box.spacing = unit(0, units = "npc"),
          legend.margin	=  margin(r = .1, l = .1, unit = "npc")) +
    NULL +
    labs(col = "% of Overlap", size = "# of Cells") +
    scale_color_gradientn(colours = RColorBrewer::brewer.pal(11, "Spectral")) +
    guides(size = guide_legend(title.position = "top", fill = "grey"),
           col = guide_colourbar(title.position = "top",
                                 barwidth = unit(.2, "npc")))
  p <- p +
    gganimate::transition_states(step,
                                 transition_length = 0,
                                 state_length =  state_length / table(Freqs$step)[1]) +
    ggtitle(paste0('Step {closest_state} of ', max(Freqs$step)))
  return(p)
}
