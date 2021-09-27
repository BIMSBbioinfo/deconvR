library(deconvR)

test_that("deconvolute nnls", {
    ref_atlas <- readRDS(system.file("reference_atlas_nodup.RDS", package = "deconvR"))
    simulation <- simulateCellMix(5)
    bulkTable <- simulation[[1]]
    proportionsTable <- simulation[[2]]

    deconvolution <- deconvolute(bulk = bulkTable, model = "nnls")
    results <- deconvolution[[1]]
    partial_rsq <- deconvolution[[2]]

    vectorlength <- NROW(tidyr::drop_na(merge(bulkTable, ref_atlas, by = "IDs"))[, -1])

    expect_output(deconvolute(bulk = bulkTable), message("DECONVOLUTION WITH NNLS"))
    expect_output(deconvolute(bulk = bulkTable), message("SUMMARY OF PARTIAL R-SQUARED VALUES FOR NNLS : "))

    expect_equal(dim(results), c(5, length(colnames(ref_atlas[, -1]))))
    expect_equal(dim(deconvolute(bulk = bulkTable[-1, ])[[1]]), c(5, length(colnames(ref_atlas[, -1]))))
    expect_equal(colnames(results), colnames(ref_atlas[, -1]))
    expect_equal(lapply(unlist(results[, ]), typeof)[[1]], print("double"),ignore_attr = TRUE)
    expect_equal(sum(results[1, ]), 1)
    expect_equal(sum(results[2, ]), 1)
    expect_equal(sum(results[3, ]), 1)
    expect_equal(sum(results[4, ]), 1)
    expect_equal(sum(results[5, ]), 1)
    expect_error(deconvolute())
    expect_error(deconvolute(bulk = bulkTable, model = "n"))
    expect_error(deconvolute(bulk = bulkTable, model = "nnls", vec = rep_len(1, vectorlength)), NA)
    expect_error(deconvolute(bulk = bulkTable, model = "nnls", vec = rep_len(1, vectorlength - 1)))
    expect_error(deconvolute(bulk = bulkTable[, -1], model = "nnls"))
    expect_error(deconvolute(model = "nnls"))
    expect_equal(dim(results), dim(proportionsTable))
    expect_equal(colnames(results), colnames(ref_atlas[, -1]))
    expect_equal(rownames(results), rownames(proportionsTable))


    expect_equal(length(partial_rsq), length(simulation[[1]]) - 1)
})

test_that("deconvolute svr", {
    ref_atlas <- readRDS(system.file("reference_atlas_nodup.RDS", package = "deconvR"))
    simulation <- simulateCellMix(5)
    bulkTable <- simulation[[1]]
    proportionsTable <- simulation[[2]]

    deconvolution <- deconvolute(bulk = bulkTable, model = "svr")
    results <- deconvolution[[1]]
    partial_rsq <- deconvolution[[2]]

    vectorlength <- NROW(tidyr::drop_na(merge(bulkTable, ref_atlas, by = "IDs"))[, -1])

    expect_output(deconvolute(bulk = bulkTable, model = "svr"), message("DECONVOLUTION WITH SVR"))
    expect_output(deconvolute(bulk = bulkTable, model = "svr"), message("SUMMARY OF PARTIAL R-SQUARED VALUES FOR SVR : "))

    ref_atlas <- readRDS(system.file("reference_atlas_nodup.RDS", package = "deconvR"))
    simulation <- simulateCellMix(5)
    bulkTable <- simulation[[1]]
    proportionsTable <- simulation[[2]]
    expect_equal(dim(results), c(5, length(colnames(ref_atlas[, -1]))))
    expect_equal(dim(deconvolute(bulk = bulkTable[-1, ], model = "svr")[[1]]), c(5, length(colnames(ref_atlas[, -1]))))
    expect_equal(colnames(results), colnames(ref_atlas[, -1]))
    expect_equal(lapply(unlist(results[, ]), typeof)[[1]], print("double"),ignore_attr = TRUE)
    expect_equal(sum(results[2, ]), 1)
    expect_equal(sum(results[3, ]), 1)
    expect_equal(sum(results[4, ]), 1)
    expect_equal(sum(results[5, ]), 1)
    expect_error(deconvolute(bulk = bulkTable[, -1], model = "svr"))
    expect_error(deconvolute(bulk = bulkTable, model = "svr", vec = rep_len(1, vectorlength)), NA)
    expect_error(deconvolute(bulk = bulkTable, model = "svr", vec = rep_len(1, vectorlength - 1)))
    expect_error(deconvolute(model = "svr"))
    expect_equal(dim(results), dim(proportionsTable))
    expect_equal(colnames(results), colnames(ref_atlas[, -1]))
    expect_equal(rownames(results), rownames(proportionsTable))

    expect_equal(length(partial_rsq), length(simulation[[1]]) - 1)
})

test_that("deconvolute qp", {
    ref_atlas <- readRDS(system.file("reference_atlas_nodup.RDS", package = "deconvR"))
    simulation <- simulateCellMix(5)
    bulkTable <- simulation[[1]]
    proportionsTable <- simulation[[2]]

    deconvolution <- deconvolute(bulk = bulkTable, model = "qp")
    results <- deconvolution[[1]]
    partial_rsq <- deconvolution[[2]]

    vectorlength <- NROW(tidyr::drop_na(merge(bulkTable, ref_atlas, by = "IDs"))[, -1])

    expect_output(deconvolute(bulk = bulkTable, model = "qp"), message("DECONVOLUTION WITH QP"))
    expect_output(deconvolute(bulk = bulkTable, model = "qp"), message("SUMMARY OF PARTIAL R-SQUARED VALUES FOR QP : "))

    expect_equal(dim(results), c(5, length(colnames(ref_atlas[, -1]))))
    expect_equal(dim(deconvolute(bulk = bulkTable[-1, ], model = "qp")[[1]]), c(5, length(colnames(ref_atlas[, -1]))))
    expect_equal(colnames(results), colnames(ref_atlas[, -1]))
    expect_equal(lapply(unlist(results[, ]), typeof)[[1]], print("double"),ignore_attr = TRUE)
    expect_equal(sum(results[1, ]), 1)
    expect_equal(sum(results[2, ]), 1)
    expect_equal(sum(results[3, ]), 1)
    expect_equal(sum(results[4, ]), 1)
    expect_equal(sum(results[5, ]), 1)
    expect_error(deconvolute(bulk = bulkTable[, -1], model = "qp"))
    expect_error(deconvolute(bulk = bulkTable, model = "qp", vec = rep_len(1, vectorlength)), NA)
    expect_error(deconvolute(bulk = bulkTable, model = "qp", vec = rep_len(1, vectorlength - 1)))
    expect_error(deconvolute(model = "qp"))
    expect_equal(dim(results), dim(proportionsTable))
    expect_equal(colnames(results), colnames(ref_atlas[, -1]))
    expect_equal(rownames(results), rownames(proportionsTable))

    expect_equal(length(partial_rsq), length(simulation[[1]]) - 1)
})

test_that("deconvolute rlm", {
    ref_atlas <- readRDS(system.file("reference_atlas_nodup.RDS", package = "deconvR"))
    simulation <- simulateCellMix(5)
    bulkTable <- simulation[[1]]
    proportionsTable <- simulation[[2]]

    deconvolution <- deconvolute(bulk = bulkTable, model = "rlm")
    results <- deconvolution[[1]]
    partial_rsq <- deconvolution[[2]]

    vectorlength <- NROW(tidyr::drop_na(merge(bulkTable, ref_atlas, by = "IDs"))[, -1])

    expect_output(deconvolute(bulk = bulkTable, model = "rlm"), message("DECONVOLUTION WITH RLM", "\n"))
    expect_output(deconvolute(bulk = bulkTable, model = "rlm"), message("SUMMARY OF PARTIAL R-SQUARED VALUES FOR RLM : "))

    expect_equal(dim(results), c(5, length(colnames(ref_atlas[, -1]))))
    expect_equal(dim(deconvolute(bulk = bulkTable[-1, ], model = "rlm")[[1]]), c(5, length(colnames(ref_atlas[, -1]))))
    expect_equal(colnames(results), colnames(ref_atlas[, -1]))
    expect_equal(lapply(unlist(results[, ]), typeof)[[1]], print("double"),ignore_attr = TRUE)
    expect_equal(sum(results[1, ]), 1)
    expect_equal(sum(results[2, ]), 1)
    expect_equal(sum(results[3, ]), 1)
    expect_equal(sum(results[4, ]), 1)
    expect_equal(sum(results[5, ]), 1)
    expect_error(deconvolute(bulk = bulkTable[, -1], model = "rlm"))
    expect_error(deconvolute(bulk = bulkTable, model = "rlm", vec = rep_len(1, vectorlength)), NA)
    expect_error(deconvolute(bulk = bulkTable, model = "rlm", vec = rep_len(1, vectorlength - 1)))
    expect_error(deconvolute(model = "rlm"))
    expect_equal(dim(results), dim(proportionsTable))
    expect_equal(colnames(results), colnames(ref_atlas[, -1]))
    expect_equal(rownames(results), rownames(proportionsTable))

    expect_equal(length(partial_rsq), length(simulation[[1]]) - 1)
})
