library(deconvR)

test_that("deconvolute nnls", {
    data("reference_atlas")
    simulation <- simulateCellMix(5, reference = reference_atlas)
    bulkTable <- simulation[[1]]
    proportionsTable <- simulation[[2]]

    deconvolution <- deconvolute(
        bulk = bulkTable,
        model = "nnls", reference = reference_atlas
    )
    results <- deconvolution[[1]]
    partial_rsq <- deconvolution[[2]]

    vectorlength <- NROW(tidyr::drop_na(merge(bulkTable, reference_atlas,
        by = "IDs"
    ))[, -1])

    expect_output(
        deconvolute(bulk = bulkTable, reference = reference_atlas),
        message("DECONVOLUTION WITH NNLS")
    )
    expect_output(
        deconvolute(bulk = bulkTable, reference = reference_atlas),
        message("SUMMARY OF PARTIAL R-SQUARED VALUES FOR NNLS : ")
    )

    expect_equal(dim(results), c(5, length(colnames(reference_atlas[, -1]))))
    expect_equal(
        dim(deconvolute(
            bulk = bulkTable[-1, ],
            reference = reference_atlas
        )[[1]]),
        c(5, length(colnames(reference_atlas[, -1])))
    )
    expect_equal(colnames(results), colnames(reference_atlas[, -1]))
    expect_equal(lapply(unlist(results[, ]), typeof)[[1]],
        print("double"),
        ignore_attr = TRUE
    )
    expect_equal(sum(results[1, ]), 1)
    expect_equal(sum(results[2, ]), 1)
    expect_equal(sum(results[3, ]), 1)
    expect_equal(sum(results[4, ]), 1)
    expect_equal(sum(results[5, ]), 1)
    expect_error(deconvolute(reference = reference_atlas))
    expect_error(deconvolute(
        bulk = bulkTable, model = "n",
        reference = reference_atlas
    ))
    expect_error(deconvolute(
        bulk = bulkTable, model = "nnls",
        vec = rep_len(1, vectorlength), reference = reference_atlas
    ), NA)
    expect_error(deconvolute(
        bulk = bulkTable, model = "nnls",
        vec = rep_len(1, vectorlength - 1), reference = reference_atlas
    ))
    expect_error(deconvolute(
        bulk = bulkTable[, -1], model = "nnls",
        reference = reference_atlas
    ))
    expect_error(deconvolute(model = "nnls", reference = reference_atlas))
    expect_equal(dim(results), dim(proportionsTable))
    expect_equal(colnames(results), colnames(reference_atlas[, -1]))
    expect_equal(rownames(results), rownames(proportionsTable))


    expect_equal(length(partial_rsq), length(simulation[[1]]) - 1)
})

test_that("deconvolute svr", {
    data("reference_atlas")
    simulation <- simulateCellMix(5, reference = reference_atlas)
    bulkTable <- simulation[[1]]
    proportionsTable <- simulation[[2]]

    deconvolution <- deconvolute(
        bulk = bulkTable, model = "svr",
        reference = reference_atlas
    )
    results <- deconvolution[[1]]
    partial_rsq <- deconvolution[[2]]

    vectorlength <- NROW(tidyr::drop_na(merge(bulkTable, reference_atlas,
        by = "IDs"
    ))[, -1])

    expect_output(
        deconvolute(
            bulk = bulkTable, model = "svr",
            reference = reference_atlas
        ),
        message("DECONVOLUTION WITH SVR")
    )
    expect_output(
        deconvolute(
            bulk = bulkTable, model = "svr",
            reference = reference_atlas
        ),
        message("SUMMARY OF PARTIAL R-SQUARED VALUES FOR SVR : ")
    )

    simulation <- simulateCellMix(5, reference = reference_atlas)
    bulkTable <- simulation[[1]]
    proportionsTable <- simulation[[2]]
    expect_equal(dim(results), c(5, length(colnames(reference_atlas[, -1]))))
    expect_equal(
        dim(deconvolute(
            bulk = bulkTable[-1, ],
            model = "svr", reference = reference_atlas
        )[[1]]),
        c(5, length(colnames(reference_atlas[, -1])))
    )
    expect_equal(colnames(results), colnames(reference_atlas[, -1]))
    expect_equal(lapply(unlist(results[, ]), typeof)[[1]],
        print("double"),
        ignore_attr = TRUE
    )
    expect_equal(sum(results[2, ]), 1)
    expect_equal(sum(results[3, ]), 1)
    expect_equal(sum(results[4, ]), 1)
    expect_equal(sum(results[5, ]), 1)
    expect_error(deconvolute(
        bulk = bulkTable[, -1], model = "svr",
        reference = reference_atlas
    ))
    expect_error(deconvolute(
        bulk = bulkTable, model = "svr",
        vec = rep_len(1, vectorlength), reference = reference_atlas
    ), NA)
    expect_error(deconvolute(
        bulk = bulkTable, model = "svr",
        vec = rep_len(1, vectorlength - 1), reference = reference_atlas
    ))
    expect_error(deconvolute(model = "svr", reference = reference_atlas))
    expect_equal(dim(results), dim(proportionsTable))
    expect_equal(colnames(results), colnames(reference_atlas[, -1]))
    expect_equal(rownames(results), rownames(proportionsTable))

    expect_equal(length(partial_rsq), length(simulation[[1]]) - 1)
})

test_that("deconvolute qp", {
    data("reference_atlas")
    simulation <- simulateCellMix(5, reference = reference_atlas)
    bulkTable <- simulation[[1]]
    proportionsTable <- simulation[[2]]

    deconvolution <- deconvolute(
        bulk = bulkTable, model = "qp",
        reference = reference_atlas
    )
    results <- deconvolution[[1]]
    partial_rsq <- deconvolution[[2]]

    vectorlength <- NROW(tidyr::drop_na(merge(bulkTable, reference_atlas,
        by = "IDs"
    ))[, -1])

    expect_output(
        deconvolute(
            bulk = bulkTable, model = "qp",
            reference = reference_atlas
        ),
        message("DECONVOLUTION WITH QP")
    )
    expect_output(
        deconvolute(
            bulk = bulkTable, model = "qp",
            reference = reference_atlas
        ),
        message("SUMMARY OF PARTIAL R-SQUARED VALUES FOR QP : ")
    )

    expect_equal(dim(results), c(5, length(colnames(reference_atlas[, -1]))))
    expect_equal(
        dim(deconvolute(
            bulk = bulkTable[-1, ], model = "qp",
            reference = reference_atlas
        )[[1]]),
        c(5, length(colnames(reference_atlas[, -1])))
    )
    expect_equal(colnames(results), colnames(reference_atlas[, -1]))
    expect_equal(lapply(unlist(results[, ]), typeof)[[1]], print("double"),
        ignore_attr = TRUE
    )
    expect_equal(sum(results[1, ]), 1)
    expect_equal(sum(results[2, ]), 1)
    expect_equal(sum(results[3, ]), 1)
    expect_equal(sum(results[4, ]), 1)
    expect_equal(sum(results[5, ]), 1)
    expect_error(deconvolute(
        bulk = bulkTable[, -1], model = "qp",
        reference = reference_atlas
    ))
    expect_error(deconvolute(
        bulk = bulkTable, model = "qp",
        vec = rep_len(1, vectorlength), reference = reference_atlas
    ), NA)
    expect_error(deconvolute(
        bulk = bulkTable, model = "qp",
        vec = rep_len(1, vectorlength - 1), reference = reference_atlas
    ))
    expect_error(deconvolute(model = "qp", reference = reference_atlas))
    expect_equal(dim(results), dim(proportionsTable))
    expect_equal(colnames(results), colnames(reference_atlas[, -1]))
    expect_equal(rownames(results), rownames(proportionsTable))

    expect_equal(length(partial_rsq), length(simulation[[1]]) - 1)
})

test_that("deconvolute rlm", {
    data("reference_atlas")
    simulation <- simulateCellMix(5, reference = reference_atlas)
    bulkTable <- simulation[[1]]
    proportionsTable <- simulation[[2]]

    deconvolution <- deconvolute(
        bulk = bulkTable, model = "rlm",
        reference = reference_atlas
    )
    results <- deconvolution[[1]]
    partial_rsq <- deconvolution[[2]]

    vectorlength <- NROW(tidyr::drop_na(merge(bulkTable, reference_atlas,
        by = "IDs"
    ))[, -1])

    expect_output(
        deconvolute(
            bulk = bulkTable, model = "rlm",
            reference = reference_atlas
        ),
        message("DECONVOLUTION WITH RLM", "\n")
    )
    expect_output(
        deconvolute(
            bulk = bulkTable, model = "rlm",
            reference = reference_atlas
        ),
        message("SUMMARY OF PARTIAL R-SQUARED VALUES FOR RLM : ")
    )

    expect_equal(dim(results), c(5, length(colnames(reference_atlas[, -1]))))
    expect_equal(
        dim(deconvolute(
            bulk = bulkTable[-1, ], model = "rlm",
            reference = reference_atlas
        )[[1]]),
        c(5, length(colnames(reference_atlas[, -1])))
    )
    expect_equal(colnames(results), colnames(reference_atlas[, -1]))
    expect_equal(lapply(unlist(results[, ]), typeof)[[1]], print("double"),
        ignore_attr = TRUE
    )
    expect_equal(sum(results[1, ]), 1)
    expect_equal(sum(results[2, ]), 1)
    expect_equal(sum(results[3, ]), 1)
    expect_equal(sum(results[4, ]), 1)
    expect_equal(sum(results[5, ]), 1)
    expect_error(deconvolute(
        bulk = bulkTable[, -1], model = "rlm",
        reference = reference_atlas
    ))
    expect_error(deconvolute(
        bulk = bulkTable, model = "rlm",
        vec = rep_len(1, vectorlength), reference = reference_atlas
    ), NA)
    expect_error(deconvolute(
        bulk = bulkTable, model = "rlm",
        vec = rep_len(1, vectorlength - 1), reference = reference_atlas
    ))
    expect_error(deconvolute(model = "rlm", reference = reference_atlas))
    expect_equal(dim(results), dim(proportionsTable))
    expect_equal(colnames(results), colnames(reference_atlas[, -1]))
    expect_equal(rownames(results), rownames(proportionsTable))

    expect_equal(length(partial_rsq), length(simulation[[1]]) - 1)
})
