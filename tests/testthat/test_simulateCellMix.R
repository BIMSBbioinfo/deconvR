library(deconvR)

test_that("simulateCellMix", {
    data("reference_atlas")
    expect_equal(
        dim(simulateCellMix(1, reference = reference_atlas)[[1]]),
        c(dim(reference_atlas[, -1])[1], 1 + 1)
    )
    expect_equal(
        simulateCellMix(1, reference = reference_atlas)[[1]][, 1],
        reference_atlas[, 1]
    )
    expect_equal(
        dim(simulateCellMix(
            numberOfSamples = 50,
            reference = reference_atlas
        )[[1]]),
        c(dim(reference_atlas[, -1])[1], 50 + 1)
    )
    expect_equal(
        simulateCellMix(
            numberOfSamples = 50,
            reference = reference_atlas
        )[[1]][, 1],
        reference_atlas[, 1]
    )
    expect_equal(sum(simulateCellMix(5,
        reference = reference_atlas
    )[[2]][1, ]), 1)
    expect_equal(sum(simulateCellMix(5,
        reference = reference_atlas
    )[[2]][2, ]), 1)
    expect_equal(sum(simulateCellMix(5,
        reference = reference_atlas
    )[[2]][3, ]), 1)
    expect_equal(sum(simulateCellMix(5,
        reference = reference_atlas
    )[[2]][4, ]), 1)
    expect_equal(sum(simulateCellMix(5,
        reference = reference_atlas
    )[[2]][5, ]), 1)
    expect_equal(sum(simulateCellMix(10,
        reference = reference_atlas
    )[[2]]), 10)
    expect_equal(
        rownames(simulateCellMix(5, reference = reference_atlas)[[2]]),
        colnames(simulateCellMix(5, reference = reference_atlas)[[1]])[-1]
    )
    expect_error(simulateCellMix(0, reference = reference_atlas)[[1]])
    expect_error(simulateCellMix(-1, reference = reference_atlas)[[1]])

    expect_equal(length(simulateCellMix(1, reference = reference_atlas, c(
        0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    ))), 2)
    expect_equal(dim(simulateCellMix(1, reference = reference_atlas, c(
        0, 0, 0, 0, 0.3, 0, 0, 0, 0.5, 0, 0,
        0, 0, 0, 0.1, 0, 0, 0, 0, 0.05, 0.05,
        0, 0, 0, 0
    ))[[1]]), c(
        nrow(reference_atlas),
        1 + 1
    ))
    expect_equal(
        dim(simulateCellMix(2, reference = reference_atlas, data.frame(c(
            0, 0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0
        ), c(
            0,
            0, 0, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0
        )))[[2]]),
        c(2, ncol(reference_atlas) - 1)
    )
    expect_equal(
        colnames(simulateCellMix(1, reference = reference_atlas, data.frame(c(
            0.1, 0.7, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0,
            0.01, 0, 0.09, 0, 0,
            0.1, 0, 0, 0, 0, 0, 0
        )))[[2]]),
        colnames(reference_atlas[, -1])
    )
    expect_error(simulateCellMix(10, reference = reference_atlas, c(
        0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    )))
    expect_error(simulateCellMix(2, reference = reference_atlas, data.frame(c(
        0.1, 0.7, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0.01, 0, 0.09, 0,
        0, 0.1, 0, 0, 0, 0, 0, 0
    ))))
    expect_error(simulateCellMix(1, reference = reference_atlas, data.frame(c(
        0, 0, 0, 0, 0.3, 0, 0, 0, 0.5,
        0, 0, 0, 0, 0, 0.1, 0, 0, 0, 0,
        0.05, 0.05, 0, 0, 0, 0
    ), c(
        0,
        0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0
    ))))
})
