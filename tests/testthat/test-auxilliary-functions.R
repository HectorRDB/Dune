context("Dune auxilliary functions")
library(Dune)
library(dplyr)

test_that("functionTracking and ARIImp are coherent", {
  data("clusMat", package = "Dune")
  merger <- Dune(clusMat)
  f <- function(clusMat) {
    ARI <- ARIs(clusMat)
    return(mean(ARI[upper.tri(ARI)]))
  }
  expect_equal(ARIImp(merger), functionTracking(merger, f))
})

test_that("clusterConversion and intermediateMat are coherent with Dune", {
  data("clusMat", package = "Dune")
  merger <- Dune(clusMat)
  df <- intermediateMat(merger, p = 1)
  df <- df[order(as.numeric(rownames(df))), ]
  expect_equal(df, merger$currentMat)
  df <- intermediateMat(merger, p = 0)
  df <- df[order(as.numeric(rownames(df))), ]
  expect_equal(df, merger$initialMat)
  df <- intermediateMat(merger, n_steps = 0)
  df <- df[order(as.numeric(rownames(df))), ]
  expect_equal(df, merger$initialMat)

})

test_that("Dune output", {
  data("clusMat", package = "Dune")
  expect_equal(sum(ARIs(Dune(clusMat)$currentMat) == 0), 0)
  expect_equal(nrow(Dune(clusMat)$merges), 7)
  expect_equal(unique(
    unlist(lapply(Dune(clusMat)$currentMat, n_distinct))),
    10)
})

test_that("Plots are returning ggplot objects", {
  data("clusMat", package = "Dune")
  merger <- Dune(clusMat)
  expect_is(plotPrePost(merger), "gg")
  expect_is(plotARIs(merger$initialMat), "gg")
  expect_is(ARItrend(merger), "gg")
})
