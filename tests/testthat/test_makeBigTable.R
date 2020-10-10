context("makeBigTable")
library(deconvR)

test_that("makeBigTable", {
  ref_atlas = read.csv(system.file("reference_atlas_nodup.csv", package = "deconvR"))
  expect_equal(dim(makeBigTable(1)), c(dim(ref_atlas[,-1])[1], 1+1))
  expect_equal(makeBigTable(1)[,1], ref_atlas[,1])
  expect_equal(dim(makeBigTable(numberOfSamples = 50, reference = ref_atlas)), c(dim(ref_atlas[,-1])[1], 50+1))
  expect_equal(makeBigTable(numberOfSamples = 50, reference = ref_atlas)[,1], ref_atlas[,1])
  expect_error(makeBigTable(0))
  expect_error(makeBigTable(-1))
})
