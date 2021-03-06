% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plotARIs}
\alias{plotARIs}
\title{Plot an heatmap of the ARI matrix}
\usage{
plotARIs(clusMat, unclustered = NULL, values = TRUE, numericalLabels = FALSE)
}
\arguments{
\item{clusMat}{The clustering matrix with a row per cell and a column per
clustering label type}

\item{unclustered}{The value assigned to unclustered cells. Default to \code{NULL}}

\item{values}{Whether to also display the ARI values. Default to TRUE.}

\item{numericalLabels}{Whether labels are numerical values. Default to FALSE.}
}
\value{
a \code{\link{ggplot}} object
}
\description{
We can compute the ARI between pairs of cluster labels. This function plots
a matrix where a cell is the adjusted Rand Index between cluster label of
row i and cluster label of column j.
}
\examples{
data("clusMat", package = "Dune")
merger <- Dune(clusMat = clusMat)
plotARIs(merger$initialMat)
plotARIs(merger$currentMat)
}
