context("Dune plotting functions do return plots")
library(Dune)
library(dplyr)

test_that("All static plotting functions return a ggplot object", {
  data("clusMat", package = "Dune")
  data("nuclei", package = "Dune")
  merger <- Dune(clusMat)
  # plotARIs 1
  p <- plotARIs(merger$initialMat)
  expect_is(p, "gg")
  expect_error(print(p), NA)
  # plotARIs 2
  p <- plotARIs(merger$initialMat, unclustered = 1)
  expect_is(p, "gg")
  expect_error(print(p), NA)
  # plotARIs 3
  mat <- merger$initialMat
  colnames(mat) <- seq_len(ncol(mat))
  p <- plotARIs(mat, numericalLabels = TRUE)
  expect_is(p, "gg")
  expect_error(print(p), NA)
  # plotPrePost
  p <- plotPrePost(merger)
  expect_is(p, "gg")
  expect_error(print(p), NA)
  # ARItrend
  p <- ARItrend(merger)
  expect_is(p, "gg")
  expect_error(print(p), NA)
  # ConfusionPlot
  p <- ConfusionPlot(nuclei[, c("SC3", "Monocle")])
  expect_is(p, "gg")
  expect_error(print(p), NA)
  # Error
  expect_error(ConfusionPlot(nuclei$SC3[1:20], nuclei$Monocle[1:10]))
})

test_that("All dynamic plotting functions return a gganim object", {
  data("clusMat", package = "Dune")
  data("nuclei", package = "Dune")
  merger <- Dune(clusMat)
  # ARI Evolution
  p <- ARIEvolution(merger)
  expect_is(ARIEvolution(merger), "gganim")
  expect_error(print(p), NA)
  # Confusion Evolution
  p <- ConfusionEvolution(merger, x = "A", y = "B")
  expect_is(p, "gganim")
  expect_error(print(p), NA)
  # Specifics
  colnames(clusMat) <- 1:5
  merger <- Dune(clusMat)
  p <- ARIEvolution(merger, numericalLabels = TRUE)
  expect_is(ARIEvolution(merger), "gganim")
  expect_error(print(p), NA)
})
