utils::globalVariables(c("new"))

#' @title ARI Matrix
#' @param clusMat The clustering matrix with a row per cell and a column per
#' clustering label type
#' @param unclustered The value assigned to unclustered cells. Default to \code{NULL}
#' @return The ARI matrix
#' @details In the ARI matrix where each cell **i,j** is the adjusted Rand Index
#' between columns **i** and **j** of the original \code{clusMat}.
#' If \code{unclustered} is not NULL, the cells which have been assigned to the
#' \code{unclustered} cluster will not be counted towards computing the ARI.
#' @importFrom aricode ARI
#' @examples
#' data("clusMat", package = "Dune")
#' ARIs(clusMat)
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
      aricode::ARI(x_unc, y_unc)
    })
  })
  if (is.null(colnames(clusMat))) {
    rownames(ARI) <- colnames(ARI) <- seq_len(ncol(clusMat))
  } else {
    rownames(ARI) <- colnames(ARI) <- colnames(clusMat)
  }

  return(ARI)
}

#' @title NMI Matrix
#' @param clusMat The clustering matrix with a row per cell and a column per
#' clustering label type
#' @param unclustered The value assigned to unclustered cells. Default to \code{NULL}
#' @return The NMI matrix
#' @details In the NMI matrix where each cell **i,j** is the normalized mutual 
#' information between columns **i** and **j** of the original \code{clusMat}.
#' If \code{unclustered} is not NULL, the cells which have been assigned to the
#' \code{unclustered} cluster will not be counted towards computing the NMI.
#' @importFrom aricode NMI
#' @examples
#' data("clusMat", package = "Dune")
#' NMIs(clusMat)
#' @export
NMIs <- function(clusMat, unclustered = NULL) {
  NMI <- apply(clusMat, 2, function(x) {
    apply(clusMat, 2, function(y) {
      if (!is.null(unclustered)) {
        x_unc <- x[x != unclustered & y != unclustered]
        y_unc <- y[x != unclustered & y != unclustered]
      } else {
        x_unc <- x
        y_unc <- y
      }
      aricode::NMI(x_unc, y_unc, variant = "sum")
    })
  })
  if (is.null(colnames(clusMat))) {
    rownames(NMI) <- colnames(NMI) <- seq_len(ncol(clusMat))
  } else {
    rownames(NMI) <- colnames(NMI) <- colnames(clusMat)
  }
  
  return(NMI)
}

#' @title When to Stop
#' @param merger the result from having run \code{\link{Dune}}
#'  on the dataset
#' @param p A value between 0 and 1. We stop when the mean ARI has improved by p
#' of the final total improvement. Default to 1 (i.e running the full merging).
#' @return An integer giving the step where to stop.
#' @details The \code{\link{Dune}} process improves the mean ARI. This return
#' the first merging step after which the mean ARI has been improved by p of the
#' total. Setting p = 1 just return the number of merges.
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

#' @title adjustedRandIndex
#' @param tab The confusion matrix
#' @return The ARI
.adjustedRandIndex <- function(tab) {
  if (all(dim(tab) == c(1, 1))) return(1)
  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2)) - a
  c <- sum(choose(colSums(tab), 2)) - a
  d <- choose(sum(tab), 2) - a - b - c
  ARI <- (a - (a + b) * (a + c)/(a + b + c + d)) /
    ((a + b + a + c) / 2 - (a + b) * (a + c) / (a + b + c + d))
  return(ARI)
}

.NormalizedMutualInformation <- function(tab){
  if (all(dim(tab) == c(1, 1))) return(1)
  freqs_x <- rowSums(tab) 
  freqs_x <- freqs_x / sum(freqs_x)
  H_x <- -sum(ifelse(freqs_x > 0, freqs_x * log(freqs_x), 0))
  freqs_y <- colSums(tab)
  freqs_y <- freqs_y / sum(freqs_y)
  H_y <- -sum(ifelse(freqs_y > 0, freqs_y * log(freqs_y), 0))
  I_xy <- base::apply(tab, 1, function(freqs) {
    freqs <- freqs / sum(freqs)
    return(-sum(ifelse(freqs > 0, freqs * log(freqs), 0)))
  })
  I_xy <- H_y - sum(freqs_x * I_xy)
  return(2 * I_xy / (H_x + H_y))
}

.confusionMatrices <- function(pairPartitions, currentMat, metric) {
  confMats <- lapply(
    as.data.frame(pairPartitions, stringsAsFactors = FALSE),
    function(pair){
      C1 <- pair[1]
      C2 <- pair[2]
      confusionMatrix <- table(currentMat[, C1], currentMat[, C2])
      if (metric == "ARI") {
        Metric <- .adjustedRandIndex(confusionMatrix)
      } else if (metric == "NMI") {
        Metric <- .NormalizedMutualInformation(confusionMatrix)
      }
      
      return(list("C1" = C1, "C2" = C2, "confusionMatrix" = confusionMatrix,
                  "Metric" = Metric))
    }
  )
  return(confMats)
}

.localARI <- function(pair, confMats, C) {
  m1 <- as.character(min(pair))
  m2 <- as.character(max(pair))
  # We look at every confusion matrix
  localARIs <- lapply(confMats, function(confMat){
    # If the confusion matrix is one between C and another partition
    if (confMat$C1 == C) {
      confusionMatrix <- confMat$confusionMatrix
      # We merge the row
      confusionMatrix[m1, ] <- confusionMatrix[m1, ] +
        confusionMatrix[m2, ]
      confusionMatrix <- confusionMatrix[
        -which(rownames(confusionMatrix) == m2), ]
      # And we update the ARI
      return(.adjustedRandIndex(confusionMatrix) - confMat$Metric)
    } else{
      # If the confusion matrix is one between another partition and C
      if (confMat$C2 == C)  {
        confusionMatrix <- confMat$confusionMatrix
        # We merge the columns
        confusionMatrix[, m1] <- confusionMatrix[, m1] +
          confusionMatrix[, m2]
        confusionMatrix <- confusionMatrix[,
                                        -which(colnames(confusionMatrix) == m2)]
        # And we update the ARI
        return(.adjustedRandIndex(confusionMatrix) - confMat$Metric)
      } else {
        # Otherwise, nothing happens in that confusion matrix
        return(0)
      }
    }
  })
  localARIs <- unlist(localARIs)
  return(localARIs)
}

.localNMI <- function(pair, confMats, C) {
  m1 <- as.character(min(pair))
  m2 <- as.character(max(pair))
  # We look at every confusion matrix
  localNMIs <- lapply(confMats, function(confMat){
    # If the confusion matrix is one between C and another partition
    if (confMat$C1 == C) {
      confusionMatrix <- confMat$confusionMatrix
      # We merge the row
      confusionMatrix[m1, ] <- confusionMatrix[m1, ] +
        confusionMatrix[m2, ]
      confusionMatrix <- confusionMatrix[
        -which(rownames(confusionMatrix) == m2), ]
      # And we update the NMI
      return(.NormalizedMutualInformation(confusionMatrix) - confMat$Metric)
    } else{
      # If the confusion matrix is one between another partition and C
      if (confMat$C2 == C)  {
        confusionMatrix <- confMat$confusionMatrix
        # We merge the columns
        confusionMatrix[, m1] <- confusionMatrix[, m1] +
          confusionMatrix[, m2]
        confusionMatrix <- confusionMatrix[,
                                           -which(colnames(confusionMatrix) == m2)]
        # And we update the NMI
        return(.NormalizedMutualInformation(confusionMatrix) - confMat$Metric)
      } else {
        # Otherwise, nothing happens in that confusion matrix
        return(0)
      }
    }
  })
  localNMIs <- unlist(localNMIs)
  return(localNMIs)
}

.updatedConfMats <- function(pair, confMats, clusLabel, metric) {
  m1 <- as.character(min(pair))
  m2 <- as.character(max(pair))
  updatedConfMats <- lapply(confMats, function(confMat){
    updatedConfMat <- confMat
    if (confMat$C1 == clusLabel) {
      confusionMatrix <- confMat$confusionMatrix
      confusionMatrix[m1, ] <- confusionMatrix[m1, ] +
        confusionMatrix[m2, ]
      confusionMatrix <- confusionMatrix[
        -which(rownames(confusionMatrix) == m2), ]
      # If we have merged all, make sure it stays in matrix format
      if (!is.matrix(confusionMatrix)) {
        confusionMatrix <- matrix(confusionMatrix, nrow = 1)
        rownames(confusionMatrix) <- m1
        colnames(confusionMatrix) <- colnames(confMat$confusionMatrix)
      }
      updatedConfMat$confusionMatrix <- confusionMatrix
      if (metric == "ARI") {
        updatedConfMat$Metric <- .adjustedRandIndex(confusionMatrix)  
      } else if (metric == "NMI") {
        updatedConfMat$Metric <- .NormalizedMutualInformation(confusionMatrix)  
      }
      
    } else{
      if (confMat$C2 == clusLabel) {
        confusionMatrix <- confMat$confusionMatrix
        confusionMatrix[, m1] <- confusionMatrix[, m1] +
          confusionMatrix[, m2]
        confusionMatrix <- confusionMatrix[,
                                           -which(colnames(confusionMatrix) == m2)]
        # If we have merged all, make sure it stays in matrix format
        if (!is.matrix(confusionMatrix)) {
          confusionMatrix <- matrix(confusionMatrix, ncol = 1)
          colnames(confusionMatrix) <- m1
          rownames(confusionMatrix) <- rownames(confMat$confusionMatrix)
        }
        updatedConfMat$confusionMatrix <- confusionMatrix
        if (metric == "ARI") {
          updatedConfMat$Metric <- .adjustedRandIndex(confusionMatrix)  
        } else if (metric == "NMI") {
          updatedConfMat$Metric <- .NormalizedMutualInformation(confusionMatrix)  
        }
      }}
    return(updatedConfMat)
  })
  return(updatedConfMats)
}

