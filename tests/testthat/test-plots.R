context("Dune plotting functions do return plots")
library(Dune)
library(dplyr)

test_that("All static plotting functions return a ggplot object", {
  data("clusMat", package = "Dune")
  data("nuclei", package = "Dune")
  merger <- Dune(clusMat)
  expect_is(plotARIs(merger$initialMat), "gg")
  expect_is(plotARIs(merger$initialMat,
                     unclustered = 1), "gg")
  expect_is(plotPrePost(merger), "gg")
  expect_is(ARItrend(merger), "gg")
  expect_is(ConfusionPlot(nuclei[, c("SC3", "Monocle")]), "gg")
})

test_that("All dynamic plotting functions return a gganim object", {
  data("clusMat", package = "Dune")
  data("nuclei", package = "Dune")
  merger <- Dune(clusMat)
  expect_is(ARIEvolution(merger), "gganim")
  expect_is(ConfusionEvolution(merger, x = "A", y = "B"), "gganim")
})
