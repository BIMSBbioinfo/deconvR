library(deconvR)
library(data.table)

test_that("findSignatures", {
  data("HumanCellTypeMethAtlas")
  # first a simple example.. only one cell type and 1 sample, no reference
  exampleSamples <- simulateCellMix(1,
    reference = HumanCellTypeMethAtlas
  )[[1]]
  exampleMeta <- data.table(
    "Experiment_accession" = "one",
    "Biosample_term_name" = "example_cell_type"
  )
  colnames(exampleSamples)[-1] <- c("one")

  results <- findSignatures(
    samples = exampleSamples,
    sampleMeta = exampleMeta,
    IDs = "IDs"
  )

  expect_equal(NCOL(results), 2)
  expect_equal(colnames(results), c("IDs", "example_cell_type"))
  expect_lte(NROW(results), NROW(exampleSamples))
  expect_equal(typeof(unlist(results[, -1])), "double")
  expect_lte(max(unlist(results[, -1])), 1)
  expect_gte(min(unlist(results[, -1])), 0)

  # new example.. three cell types and 25 samples, still no reference
  data("HumanCellTypeMethAtlas")
  exampleSamples <- HumanCellTypeMethAtlas
  exampleMeta <- data.table(
    "Experiment_accession" = as.character(1:25),
    "Biosample_term_name" = c(
      rep("this", 2), rep("that", 3), "this",
      rep("else", 14), rep("that", 5)
    )
  )
  exampleMeta <- exampleMeta[sample(nrow(exampleMeta)), ]
  colnames(exampleSamples)[-1] <- as.character(1:25)
  results <- findSignatures(
    samples = exampleSamples,
    sampleMeta = exampleMeta,
    IDs = "IDs"
  )

  expect_equal(NCOL(results), 4)
  expect_setequal(colnames(results), c("IDs", "this", "that", "else"))
  expect_equal(colnames(results)[1], "IDs")
  expect_lte(NROW(results), NROW(exampleSamples))
  expect_equal(typeof(unlist(results[, -1])), "double")
  expect_lte(max(unlist(results[, -1])), 1)
  expect_gte(min(unlist(results[, -1])), 0)

  results_cutoff <- findSignatures(
    samples = exampleSamples,
    sampleMeta = exampleMeta,
    variation_cutoff = 0.01,
    IDs = "IDs"
  )
  expect_gt(NROW(results), NROW(results_cutoff))


  # now a more complex example, with 10 cell types, 100 samples,
  # and a reference atlas of 25 cell types
  data("HumanCellTypeMethAtlas")
  exampleSamples <- simulateCellMix(100,
    reference =
      HumanCellTypeMethAtlas[c(1:1000, 2000:5000), ]
  )[[1]]
  HumanCellTypeMethAtlas <-
    HumanCellTypeMethAtlas[sample(nrow(HumanCellTypeMethAtlas)), ]
  exampleMeta <- data.table(
    "Experiment_accession" = c(
      colnames(exampleSamples)[-1],
      colnames(HumanCellTypeMethAtlas)[-1],
      as.character(1:200)
    ),
    "Biosample_term_name" = c(
      rep("this", 2), rep("that", 13), rep("else", 25),
      rep(colnames(HumanCellTypeMethAtlas)[2], 1),
      rep("this_v2", 31), rep("that_v2", 17),
      rep("else_v2", 90), rep("other", 146)
    )
  )
  exampleMeta <- exampleMeta[sample(nrow(exampleMeta)), ]
  results <- findSignatures(
    samples = exampleSamples,
    sampleMeta = exampleMeta,
    atlas = HumanCellTypeMethAtlas,
    IDs = "IDs"
  )

  expect_equal(NCOL(results), 8)
  expect_setequal(colnames(results), c(
    "IDs", "this", "that", "else",
    "this_v2", "that_v2", "else_v2",
    colnames(HumanCellTypeMethAtlas)[2]
  ))
  expect_equal(colnames(results)[1], "IDs")
  expect_lte(NROW(results), NROW(exampleSamples))
  expect_equal(typeof(unlist(results[, -1])), "double")
  expect_lte(max(unlist(results[, -1])), 1)
  expect_gte(min(unlist(results[, -1])), 0)

  ## construct tissue specific methylation signatures
  data("HumanCellTypeMethAtlas")
  exampleSamples <- simulateCellMix(1,
    reference = HumanCellTypeMethAtlas
  )$simulated
  exampleMeta <- data.table(
    "Experiment_accession" = "example_sample",
    "Biosample_term_name" = "example_cell_type"
  )
  colnames(exampleSamples) <- c("CpGs", "example_sample")
  colnames(HumanCellTypeMethAtlas)[1] <- c("CpGs")
  signatures <- findSignatures(
    samples = exampleSamples,
    sampleMeta = exampleMeta,
    atlas = HumanCellTypeMethAtlas,
    IDs = "CpGs", K = 100, tissueSpecCpGs = TRUE
  )
  expect_equal(names(signatures)[1], "example_cell_type")
  expect_lte(NROW(signatures), NCOL(HumanCellTypeMethAtlas))
  expect_equal(typeof(unlist(signatures)), "double")
  expect_lte(max(unlist(signatures)), 1)
  expect_gte(min(unlist(signatures)), 0)

  signatures <- findSignatures(
    samples = exampleSamples,
    sampleMeta = exampleMeta,
    atlas = HumanCellTypeMethAtlas,
    IDs = "CpGs", tissueSpecDMPs = TRUE
  )
  expect_equal(names(signatures)[1], "example_sample")
  expect_lte(NROW(signatures), NCOL(HumanCellTypeMethAtlas))
  expect_equal(typeof(unlist(signatures)), "double")
})
