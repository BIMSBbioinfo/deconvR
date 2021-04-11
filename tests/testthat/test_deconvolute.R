context("deconvolution")
library(deconvR)

test_that("deconvolute nnls", {
  ref_atlas = readRDS(system.file("reference_atlas_nodup.RDS", package = "deconvR"))
  simulation = simulateCellMix(5)
  bulkTable = simulation[[1]]
  proportionsTable = simulation[[2]]

  results = deconvolute(bulk=bulkTable, model="nnls")

  expect_output(deconvolute(bulk=bulkTable), "DECONVOLUTION WITH NNLS")
  expect_output(deconvolute(bulk=bulkTable), "SUMMARY OF PARTIAL R-SQUARED VALUES FOR NNLS : ")

  expect_equal(dim(results), c(5,length(colnames(ref_atlas[,-1]))))
  expect_equal(dim(deconvolute(bulk=bulkTable[-1,])), c(5,length(colnames(ref_atlas[,-1]))))
  expect_equal(colnames(results), colnames(ref_atlas[,-1]))
  expect_equivalent(sapply(results[,], typeof), rep_len("double", length(colnames(ref_atlas[,-1]))))
  expect_equal(sum(results[1,]), 1)
  expect_equal(sum(results[2,]), 1)
  expect_equal(sum(results[3,]), 1)
  expect_equal(sum(results[4,]), 1)
  expect_equal(sum(results[5,]), 1)
  expect_error(deconvolute())
  expect_error(deconvolute(bulk=bulkTable, model="n"))
  expect_error(deconvolute(bulk=bulkTable[,-1], model="nnls"))
  expect_error(deconvolute(model="nnls"))
  expect_equal(dim(results), dim(proportionsTable))
  expect_equal(colnames(results), colnames(ref_atlas[,-1]))
  expect_equal(rownames(deconvolute(bulk=bulkTable, model ="nnls")), rownames(proportionsTable))
})

test_that("deconvolute svr", {
  ref_atlas = readRDS(system.file("reference_atlas_nodup.RDS", package = "deconvR"))
  simulation = simulateCellMix(5)
  bulkTable = simulation[[1]]
  proportionsTable = simulation[[2]]

  results = deconvolute(bulk=bulkTable, model="svr")

  expect_output(deconvolute(bulk=bulkTable, model="svr"), "DECONVOLUTION WITH SVR")
  expect_output(deconvolute(bulk=bulkTable, model="svr"), "SUMMARY OF PARTIAL R-SQUARED VALUES FOR SVR : ")

  ref_atlas = readRDS(system.file("reference_atlas_nodup.RDS", package = "deconvR"))
  simulation = simulateCellMix(5)
  bulkTable = simulation[[1]]
  proportionsTable = simulation[[2]]
  expect_equal(dim(results), c(5,length(colnames(ref_atlas[,-1]))))
  expect_equal(dim(deconvolute(bulk=bulkTable[-1,], model="svr")), c(5,length(colnames(ref_atlas[,-1]))))
  expect_equal(colnames(results), colnames(ref_atlas[,-1]))
  expect_equivalent(sapply(results[,], typeof), rep_len("double", length(colnames(ref_atlas[,-1]))))
  expect_equal(sum(results[2,]), 1)
  expect_equal(sum(results[3,]), 1)
  expect_equal(sum(results[4,]), 1)
  expect_equal(sum(results[5,]), 1)
  expect_error(deconvolute(bulk=bulkTable[,-1], model="svr"))
  expect_error(deconvolute(model="svr"))
  expect_equal(dim(results), dim(proportionsTable))
  expect_equal(colnames(results), colnames(ref_atlas[,-1]))
  expect_equal(rownames(deconvolute(bulk=bulkTable, model ="svr")), rownames(proportionsTable))
})

test_that("deconvolute qp", {
  ref_atlas = readRDS(system.file("reference_atlas_nodup.RDS", package = "deconvR"))
  simulation = simulateCellMix(5)
  bulkTable = simulation[[1]]
  proportionsTable = simulation[[2]]

  results = deconvolute(bulk=bulkTable, model="qp")

  expect_output(deconvolute(bulk=bulkTable, model="qp"), "DECONVOLUTION WITH QP")
  expect_output(deconvolute(bulk=bulkTable, model="qp"), "SUMMARY OF PARTIAL R-SQUARED VALUES FOR QP : ")

  expect_equal(dim(results), c(5,length(colnames(ref_atlas[,-1]))))
  expect_equal(dim(deconvolute(bulk=bulkTable[-1,], model="qp")), c(5,length(colnames(ref_atlas[,-1]))))
  expect_equal(colnames(results), colnames(ref_atlas[,-1]))
  expect_equivalent(sapply(results[,], typeof), rep_len("double", length(colnames(ref_atlas[,-1]))))
  expect_equal(sum(results[1,]), 1)
  expect_equal(sum(results[2,]), 1)
  expect_equal(sum(results[3,]), 1)
  expect_equal(sum(results[4,]), 1)
  expect_equal(sum(results[5,]), 1)
  expect_error(deconvolute(bulk=bulkTable[,-1], model="qp"))
  expect_error(deconvolute(model="qp"))
  expect_equal(dim(results), dim(proportionsTable))
  expect_equal(colnames(results), colnames(ref_atlas[,-1]))
  expect_equal(rownames(deconvolute(bulk=bulkTable, model ="qp")), rownames(proportionsTable))
})

test_that("deconvolute rlm", {
  ref_atlas = readRDS(system.file("reference_atlas_nodup.RDS", package = "deconvR"))
  simulation = simulateCellMix(5)
  bulkTable = simulation[[1]]
  proportionsTable = simulation[[2]]

  results = deconvolute(bulk=bulkTable, model="rlm")

  expect_output(deconvolute(bulk=bulkTable, model="rlm"), "DECONVOLUTION WITH RLM")
  expect_output(deconvolute(bulk=bulkTable, model="rlm"), "SUMMARY OF PARTIAL R-SQUARED VALUES FOR RLM : ")

  expect_equal(dim(results), c(5,length(colnames(ref_atlas[,-1]))))
  expect_equal(dim(deconvolute(bulk=bulkTable[-1,], model="rlm")), c(5,length(colnames(ref_atlas[,-1]))))
  expect_equal(colnames(results), colnames(ref_atlas[,-1]))
  expect_equivalent(sapply(results[,], typeof), rep_len("double", length(colnames(ref_atlas[,-1]))))
  expect_equal(sum(results[1,]), 1)
  expect_equal(sum(results[2,]), 1)
  expect_equal(sum(results[3,]), 1)
  expect_equal(sum(results[4,]), 1)
  expect_equal(sum(results[5,]), 1)
  expect_error(deconvolute(bulk=bulkTable[,-1], model="rlm"))
  expect_error(deconvolute(model="rlm"))
  expect_equal(dim(results), dim(proportionsTable))
  expect_equal(colnames(results), colnames(ref_atlas[,-1]))
  expect_equal(rownames(deconvolute(bulk=bulkTable, model ="rlm")), rownames(proportionsTable))

})
