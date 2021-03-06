context("Dune's main algorithm works as intended")

library(Dune)
library(dplyr)
set.seed(107070770)

test_that("Dune returns the right type of output", {
  data("clusMat", package = "Dune")
  merger <- Dune(clusMat)
  expect_equal(length(merger), 5)
  expect_equal(names(merger),
               c("initialMat", "currentMat", "merges", "ImpMetric", "metric"))
  expect_is(merger$initialMat, "matrix")
  expect_is(merger$currentMat, "data.frame")
  expect_is(merger$merges, "data.frame")
  expect_is(merger$ImpMetric, "numeric")
  expect_output(Dune(clusMat, verbose = TRUE))
})


test_that("Dune correctly picks the best cluster for ARI", {
  for (i in 1:10) {
    clusMat <- matrix(sample(1:5, 500, replace = TRUE), ncol = 5)
    merger <- Dune(clusMat, metric = "ARI")
    if (nrow(merger$merges) <= 3) next()
    df <- intermediateMat(merger, n_steps = nrow(merger$merges) - 2)
    df <- as.matrix(df[, -1])
    init_ARI <- ARIs(df)
    for (j in 1:20) {
      col <- sample(colnames(df), 1)
      pair <- sample(t(unique(df[,col])), 2)
      m2 <- max(pair)
      m1 <- min(pair)
      df2 <- df
      df2[df2[, col] == m2, col] <- m1
      final_ARI <- ARIs(df2)
      expect_true(mean((final_ARI - init_ARI)[upper.tri(init_ARI)]) <=
                    merger$ImpMetric[nrow(merger$merges) - 1])
    }
  }
})

test_that("Dune correctly picks the best cluster for NMI", {
  for (i in 1:10) {
    clusMat <- matrix(sample(1:5, 500, replace = TRUE), ncol = 5)
    merger <- Dune(clusMat, metric = "ARI")
    if (nrow(merger$merges) <= 3) next()
    df <- intermediateMat(merger, n_steps = nrow(merger$merges) - 2)
    df <- as.matrix(df[, -1])
    init_ARI <- ARIs(df)
    for (j in 1:20) {
      col <- sample(colnames(df), 1)
      pair <- sample(t(unique(df[,col])), 2)
      m2 <- max(pair)
      m1 <- min(pair)
      df2 <- df
      df2[df2[, col] == m2, col] <- m1
      final_ARI <- ARIs(df2)
      expect_true(mean((final_ARI - init_ARI)[upper.tri(init_ARI)]) <=
                    merger$ImpMetric[nrow(merger$merges) - 1])
    }
  }
})

test_that('Dune works correctly for identical inputs',{
  data("clusMat", package = "Dune")
  expect_error(Dune(clusMat[,c(1,1)]))
  expect_output(Dune(cbind(rep(1:3, each = 20), rep(1:2, each = 30)),
                     verbose = TRUE),
                regexp = ".*We merge one of the partition entirely")
})


test_that('Dune works with unclustered values',{
  data("clusMat", package = "Dune")
  merger <- Dune(clusMat, unclustered = -2)
  merger2 <- Dune(clusMat)
  expect_identical(merger, merger2)
})
