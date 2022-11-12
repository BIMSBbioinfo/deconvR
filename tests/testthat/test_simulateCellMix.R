library(deconvR)

test_that("simulateCellMix", {
  data("HumanCellTypeMethAtlas")
  expect_equal(
    dim(simulateCellMix(1, reference = HumanCellTypeMethAtlas)$simulated),
    c(dim(HumanCellTypeMethAtlas[, -1])[1], 1 + 1)
  )
  expect_equal(
    simulateCellMix(1, reference = HumanCellTypeMethAtlas)$simulated[, 1],
    HumanCellTypeMethAtlas[, 1]
  )
  expect_equal(
    dim(simulateCellMix(
      numberOfSamples = 50,
      reference = HumanCellTypeMethAtlas
    )$simulated),
    c(dim(HumanCellTypeMethAtlas[, -1])[1], 50 + 1)
  )
  expect_equal(
    simulateCellMix(
      numberOfSamples = 50,
      reference = HumanCellTypeMethAtlas
    )$simulated[, 1],
    HumanCellTypeMethAtlas[, 1]
  )
  expect_equal(sum(simulateCellMix(5,
    reference = HumanCellTypeMethAtlas
  )[[2]][1, ]), 1)
  expect_equal(sum(simulateCellMix(5,
    reference = HumanCellTypeMethAtlas
  )[[2]][2, ]), 1)
  expect_equal(sum(simulateCellMix(5,
    reference = HumanCellTypeMethAtlas
  )[[2]][3, ]), 1)
  expect_equal(sum(simulateCellMix(5,
    reference = HumanCellTypeMethAtlas
  )[[2]][4, ]), 1)
  expect_equal(sum(simulateCellMix(5,
    reference = HumanCellTypeMethAtlas
  )[[2]][5, ]), 1)
  expect_equal(sum(simulateCellMix(10,
    reference = HumanCellTypeMethAtlas
  )[[2]]), 10)
  expect_equal(
    rownames(simulateCellMix(5,
      reference = HumanCellTypeMethAtlas
    )[[2]]),
    colnames(simulateCellMix(5,
      reference = HumanCellTypeMethAtlas
    )$simulated)[-1]
  )
  expect_error(simulateCellMix(0,
    reference = HumanCellTypeMethAtlas
  )$simulated)
  expect_error(simulateCellMix(-1,
    reference = HumanCellTypeMethAtlas
  )$simulated)

  expect_equal(length(simulateCellMix(1,
    reference = HumanCellTypeMethAtlas, c(
      0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    )
  )), 2)
  expect_equal(dim(simulateCellMix(1,
    reference = HumanCellTypeMethAtlas, c(
      0, 0, 0, 0, 0.3, 0, 0, 0, 0.5, 0, 0,
      0, 0, 0, 0.1, 0, 0, 0, 0, 0.05, 0.05,
      0, 0, 0, 0
    )
  )$simulated), c(
    nrow(HumanCellTypeMethAtlas),
    1 + 1
  ))
  expect_equal(
    dim(simulateCellMix(2,
      reference = HumanCellTypeMethAtlas, data.frame(c(
        0, 0, 0, 0, 0, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0
      ), c(
        0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0
      ))
    )[[2]]),
    c(2, ncol(HumanCellTypeMethAtlas) - 1)
  )
  expect_equal(
    colnames(simulateCellMix(1,
      reference = HumanCellTypeMethAtlas, data.frame(c(
        0.1, 0.7, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0,
        0.01, 0, 0.09, 0, 0,
        0.1, 0, 0, 0, 0, 0, 0
      ))
    )[[2]]),
    colnames(HumanCellTypeMethAtlas[, -1])
  )
  expect_error(simulateCellMix(10, reference = HumanCellTypeMethAtlas, c(
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
  )))
  expect_error(simulateCellMix(2,
    reference = HumanCellTypeMethAtlas, data.frame(c(
      0.1, 0.7, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0.01, 0, 0.09, 0,
      0, 0.1, 0, 0, 0, 0, 0, 0
    ))
  ))
  expect_error(simulateCellMix(1,
    reference = HumanCellTypeMethAtlas,
    data.frame(c(
      0, 0, 0, 0, 0.3, 0, 0, 0, 0.5,
      0, 0, 0, 0, 0, 0.1, 0, 0, 0, 0,
      0.05, 0.05, 0, 0, 0, 0
    ), c(
      0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0
    ))
  ))
})
