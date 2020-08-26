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
  expect_equal(ARIImp(merger),
               functionTracking(merger, f, n_steps = length(merger$ImpMetric)))
})

test_that("clusterConversion and intermediateMat are coherent with Dune", {
  data("clusMat", package = "Dune")
  merger <- Dune(clusMat)
  # Same at the end 
  df <- intermediateMat(merger, p = 1) %>%
    mutate(cells = as.numeric(cells)) %>%
    arrange(cells) %>%
    select(-cells) %>%
    as.data.frame()
  expect_equal(df, merger$currentMat)
  # Same at the beginning
  df <- intermediateMat(merger, p = 0) %>%
    mutate(cells = as.numeric(cells)) %>%
    arrange(cells) %>%
    select(-cells) %>%
    as.matrix()
  expect_equal(df, merger$initialMat)
  # Same at the beginning
  df <- intermediateMat(merger, n_steps = 0) %>%
    mutate(cells = as.numeric(cells)) %>%
    arrange(cells) %>%
    select(-cells) %>%
    as.matrix()
  expect_equal(df, merger$initialMat)
  # Same if rownames
  rownames(merger$initialMat) <- seq_len(nrow(clusMat))
  df <- intermediateMat(merger, p = 1) %>%
    mutate(cells = as.numeric(cells)) %>%
    arrange(cells) %>%
    select(-cells) %>%
    as.data.frame()
  expect_equal(df, merger$currentMat)
})

test_that("Dune output", {
  data("clusMat", package = "Dune")
  expect_equal(sum(ARIs(Dune(clusMat)$currentMat) == 0), 0)
  expect_equal(nrow(Dune(clusMat)$merges), 7)
  expect_equal(unique(
    unlist(lapply(Dune(clusMat)$currentMat, n_distinct))),
    10)
})

test_that("You can only stop between zero and one", {
  data("clusMat", package = "Dune")
  merger <- Dune(clusMat)
  expect_error(whenToStop(merger, -1))
  expect_error(whenToStop(merger, 2))
  expect_error(whenToStop(merger, "A"))
})


