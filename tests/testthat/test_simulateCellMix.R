context("simulateCellMix")
library(deconvR)

test_that("simulateCellMix", {
  ref_atlas = readRDS(system.file("reference_atlas_nodup.RDS", package = "deconvR"))
  expect_equal(dim(simulateCellMix(1)[[1]]), c(dim(ref_atlas[,-1])[1], 1+1))
  expect_equal(simulateCellMix(1)[[1]][,1], ref_atlas[,1])
  expect_equal(dim(simulateCellMix(numberOfSamples = 50, reference = ref_atlas)[[1]]), c(dim(ref_atlas[,-1])[1], 50+1))
  expect_equal(simulateCellMix(numberOfSamples = 50, reference = ref_atlas)[[1]][,1], ref_atlas[,1])
  expect_error(simulateCellMix(0)[[1]])
  expect_error(simulateCellMix(-1)[[1]])

  expect_equal(length(simulateCellMix(1, c(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))), 2)
  expect_equal(dim(simulateCellMix(1, c(0,0,0,0,0.3,0,0,0,0.5,0,0,0,0,0,0.1,0,0,0,0,0.05,0.05,0,0,0,0))[[1]]), c(nrow(ref_atlas),1+1))
  expect_equal(dim(simulateCellMix(2, data.frame(c(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), c(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)))[[2]]), c(2,ncol(ref_atlas)-1))
  expect_equal(colnames(simulateCellMix(1, data.frame(c(0.1,0.7,0,0,0,0,0,0,0,0,0,0,0,0.01,0,0.09,0,0,0.1,0,0,0,0,0,0)))[[2]]), colnames(ref_atlas[,-1]))
  expect_error(simulateCellMix(10, c(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)))
  expect_error(simulateCellMix(2, data.frame(c(0.1,0.7,0,0,0,0,0,0,0,0,0,0,0,0.01,0,0.09,0,0,0.1,0,0,0,0,0,0))))
  expect_error(simulateCellMix(1, data.frame(c(0,0,0,0,0.3,0,0,0,0.5,0,0,0,0,0,0.1,0,0,0,0,0.05,0.05,0,0,0,0),  c(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))))

  expect_equivalent(data.frame(t(c(0.1,0.7,0,0,0,0,0,0,0,0,0,0,0,0.01,0,0.09,0,0,0.1,0,0,0,0,0,0))), simulateCellMix(1, data.frame(c(0.1,0.7,0,0,0,0,0,0,0,0,0,0,0,0.01,0,0.09,0,0,0.1,0,0,0,0,0,0)))[[2]])
  expect_equivalent(data.frame(rbind(c(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), c(0,0,0,0,0.3,0,0,0,0.5,0,0,0,0,0,0.1,0,0,0,0,0.05,0.05,0,0,0,0))), simulateCellMix(2, data.frame(c(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),  c(0,0,0,0,0.3,0,0,0,0.5,0,0,0,0,0,0.1,0,0,0,0,0.05,0.05,0,0,0,0)))[[2]])
  expect_equivalent(data.frame(rbind(c(0.1,0.7,0,0,0,0,0,0,0,0,0,0,0,0.01,0,0.09,0,0,0.1,0,0,0,0,0,0))), simulateCellMix(1, c(0.1,0.7,0,0,0,0,0,0,0,0,0,0,0,0.01,0,0.09,0,0,0.1,0,0,0,0,0,0))[[2]])

})