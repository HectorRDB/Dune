library(Dune)
library(dplyr)

test_that("Dune returns the right type of output", {
  data("clusMat", package = "Dune")
  merger <- Dune(clusMat)
  expect_equal(length(merger), 4)
  expect_equal(names(merger), c("initialMat", "currentMat", "merges", "ImpARI"))
  expect_is(merger$initialMat, "matrix")
  expect_is(merger$currentMat, "data.frame")
  expect_is(merger$merges, "data.frame")
  expect_is(merger$ImpARI, "numeric")
})

test_that("Dune correctly picks the best cluster", {
  for (i in 1:10) {
    clusMat <- matrix(sample(1:5, 500, replace = TRUE), ncol = 5)
    merger <- Dune(clusMat)
    df <- intermediateMat(merger, n_steps = nrow(merger$merges) - 2)
    df <- as.matrix(df)
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
                    merger$ImpARI[nrow(merger$merges) - 1])
    }
  }
})
