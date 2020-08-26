context("Dune works as expected for all inputs")
library(Dune)
library(SummarizedExperiment)

test_that("Dune reacts similarly to all kind of inputs", {
  data("clusMat", package = "Dune")
  merger_matrix <- Dune(clusMat)
  merger_df <- Dune(as.data.frame(clusMat))
  sce <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = t(clusMat)),
                                                    colData = clusMat)
  merger_sExp <- Dune(clusMat = sce, cluster_columns = colnames(clusMat))
  expect_identical(merger_matrix$initialMat, as.matrix(merger_df$initialMat))
  expect_identical(merger_matrix$initialMat, as.matrix(merger_sExp$initialMat))
  expect_identical(merger_matrix[2:4], merger_df[2:4])
  expect_identical(merger_matrix[2:4], merger_sExp[2:4])
})

test_that("Dune fails when wrong column names are specified", {
  data("clusMat", package = "Dune")
  sce <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = t(clusMat)),
                                                    colData = clusMat)
  expect_error(Dune(clusMat = sce, cluster_columns = c("A", "F")))
})

test_that("Dune fails when the wrong metric is specified", {
  data("clusMat", package = "Dune")
  expect_error(Dune(clusMat = clusMat, metric = "AEIHn"))
})