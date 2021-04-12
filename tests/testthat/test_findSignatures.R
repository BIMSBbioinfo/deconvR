context("findSignatures")
library(deconvR)

test_that("findSignatures", {
  # first a simple example.. only one cell type and 1 sample, no reference
  exampleSamples = simulateCellMix(1)[[1]]
  exampleMeta  = data.table("Experiment_accession" = "one",
                            "Biosample_term_name" = "example_cell_type")
  colnames(exampleSamples)[-1] =  c("one")

  results = findSignatures(samples = exampleSamples, sampleMeta = exampleMeta)

  expect_equal(NCOL(results), 2)
  expect_equal(colnames(results), c("IDs", "example_cell_type"))
  expect_lte(NROW(results), NROW(exampleSamples))
  expect_equivalent(typeof(unlist(results[,-1])), "double")
  expect_lte(max(unlist(results[,-1])), 1)
  expect_gte(min(unlist(results[,-1])), 0)

  # new example.. three cell types and 25 samples, still no reference
  exampleSamples = readRDS(system.file("reference_atlas_nodup.RDS", package = "deconvR"))
  exampleMeta  = data.table("Experiment_accession" = as.character(1:25),
                            "Biosample_term_name" = c(rep("this",2), rep("that", 3), "this", rep("else", 14), rep("that", 5)))
  exampleMeta = exampleMeta[sample(nrow(exampleMeta)),]
  colnames(exampleSamples)[-1] =  as.character(1:25)
  results = findSignatures(samples = exampleSamples, sampleMeta = exampleMeta)

  expect_equal(NCOL(results), 4)
  expect_setequal(colnames(results), c("IDs", "this", "that", "else"))
  expect_equal(colnames(results)[1], "IDs")
  expect_lte(NROW(results), NROW(exampleSamples))
  expect_equivalent(typeof(unlist(results[,-1])), "double")
  expect_lte(max(unlist(results[,-1])), 1)
  expect_gte(min(unlist(results[,-1])), 0)

  results_cutoff = findSignatures(samples = exampleSamples, sampleMeta = exampleMeta, variation_cutoff = 0.01)
  expect_gt(NROW(results), NROW(results_cutoff))


  # now a more complex example, with 10 cell types, 100 samples, and a reference atlas of 25 cell types
  exampleReference = readRDS(system.file("reference_atlas_nodup.RDS", package = "deconvR"))
  exampleSamples = simulateCellMix(100, reference = exampleReference[c(1:1000, 2000:5000),])[[1]]
  exampleReference = exampleReference[sample(nrow(exampleReference)),]
  exampleMeta  = data.table("Experiment_accession" = c(colnames(exampleSamples)[-1],colnames(exampleReference)[-1], as.character(1:200)),
                            "Biosample_term_name" = c(rep("this",2), rep("that", 13), rep("else", 25), rep(colnames(exampleReference)[2], 1),
                                                      rep("this_v2",31), rep("that_v2", 17), rep("else_v2", 90), rep("other", 146)))
  exampleMeta = exampleMeta[sample(nrow(exampleMeta)),]
  results = findSignatures(samples = exampleSamples, sampleMeta = exampleMeta, atlas = exampleReference)

  expect_equal(NCOL(results), 8)
  expect_setequal(colnames(results), c("IDs", "this", "that", "else","this_v2", "that_v2", "else_v2", colnames(exampleReference)[2]))
  expect_equal(colnames(results)[1], "IDs")
  expect_lte(NROW(results), NROW(exampleSamples))
  expect_equivalent(typeof(unlist(results[,-1])), "double")
  expect_lte(max(unlist(results[,-1])), 1)
  expect_gte(min(unlist(results[,-1])), 0)

  })

