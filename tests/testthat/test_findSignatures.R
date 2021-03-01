context("findSignatures")
library(deconvR)

test_that("findSignatures", {
  exampleAtlas = simulateCellMix(5)[[1]]
  exampleMeta  = data.table("Experiment_accession" = c("one", "two", "three", "four", "five"), 
                            "Biosample_term_name" = "example_cell_type")
  colnames(exampleAtlas)[-1] =  c("one", "two", "three", "four", "five")
  
  
  expect_equal(NCOL(findSignatures(samples = exampleAtlas, sampleMeta = exampleMeta)), 2)
  expect_equal(colnames(findSignatures(samples = exampleAtlas, sampleMeta = exampleMeta)), c("CpGs", "example_cell_type"))
  expect_lte(NROW(findSignatures(samples = exampleAtlas, sampleMeta = exampleMeta)), NROW(exampleAtlas))
  
  
  
  exampleAtlas = readRDS(system.file("reference_atlas_nodup.RDS", package = "deconvR"))
  exampleMeta  = data.table("Experiment_accession" = as.character(1:25), 
                            "Biosample_term_name" = "example_cell_type")
  colnames(exampleAtlas)[-1] =  as.character(1:25)
  expect_equal(findSignatures(samples = exampleAtlas, sampleMeta = exampleMeta), exampleAtlas)
  
  exampleAtlas_2 = readRDS(system.file("reference_atlas_nodup.RDS", package = "deconvR"))
  exampleMeta_2  = data.table("Experiment_accession" = as.character(1:25), "Biosample_term_name" = colnames(exampleAtlas_2)[-1])
  exampleSamples_2 = exampleAtlas_2[,c(1,26)]
  exampleAtlas_2 = exampleAtlas_2[,c(1:25)]
  exampleMeta_2    = exampleMeta_2[25,]
  colnames(exampleSamples_2)[2]   = "25"
  findSignatures(samples = exampleSamples_2, sampleMeta = exampleMeta_2, atlas = exampleAtlas_2)
  
  })

