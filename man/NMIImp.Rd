% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core.R
\name{NMIImp}
\alias{NMIImp}
\title{NMI improvement}
\usage{
NMIImp(merger, unclustered = NULL)
}
\arguments{
\item{merger}{the result from having run \code{\link{Dune}}
on the dataset}

\item{unclustered}{The value assigned to unclustered cells. Default to \code{NULL}}
}
\value{
a vector with the mean NMI between methods at each merge
}
\description{
Compute the NMI improvement over the NMI merging procedure
}
\examples{
data("clusMat", package = "Dune")
merger <- Dune(clusMat = clusMat)
plot(0:nrow(merger$merges), NMIImp(merger))
}
\seealso{
NMItrend
}
