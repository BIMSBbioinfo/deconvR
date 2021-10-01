## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(tidy = FALSE,
                      cache = FALSE,
                      dev = "png",
                      message = FALSE, 
                      error = FALSE,
                      warning = FALSE)
BiocStyle::markdown()
library(knitr)
library(deconvR)
library(doParallel)
library(dplyr)

cl <- parallel::makeCluster(2)
doParallel::registerDoParallel(cl)

## ----eval=FALSE---------------------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  
#  BiocManager::install("deconvR")

## ----eval=FALSE---------------------------------------------------------------
#  simulateCellMix(
#    numberOfSamples,
#    mixingVector=NULL,
#    reference = reference_atlas)

## ----eval=FALSE---------------------------------------------------------------
#  deconvolute(
#    reference = reference_atlas,vec = NULL, bulk,
#    model= "nnls")
#  

## ----eval=FALSE---------------------------------------------------------------
#  BSmeth2Probe(
#    probe_id_locations,
#    WGBS_data,
#    cutoff = 10,
#    multipleMapping = FALSE)

## ----eval=FALSE---------------------------------------------------------------
#  findSignatures(samples,
#                 sampleMeta,
#                 atlas = NULL,
#                 variation_cutoff = NULL)

## ---- message = FALSE, output.lines=10----------------------------------------
data("WGBS_GRanges")
WGBS_GRanges

## ---- message = FALSE, output.lines=10----------------------------------------
data("probe_ids")
probe_ids

## ---- message = FALSE, output.lines=10----------------------------------------
library(deconvR) 

data("reference_atlas")
head(reference_atlas)

## ---- message = FALSE, output.lines=10----------------------------------------
samples <- simulateCellMix(3,reference = reference_atlas)[[1]]
head(samples)

## ---- message = FALSE, output.lines=10----------------------------------------
sampleMeta <- data.table("Experiment_accession" = colnames(samples)[-1],
                         "Biosample_term_name" = "new cell type")
head(sampleMeta)

## ---- output.lines=10---------------------------------------------------------
extended_matrix <- findSignatures(samples = samples, 
                                 sampleMeta = sampleMeta, 
                                 atlas = reference_atlas)
head(extended_matrix)

## ---- message = FALSE, output.lines=10----------------------------------------
data("WGBS_GRanges")
WGBS_GRanges

## ---- message = FALSE, output.lines=10----------------------------------------
head(methylKit::methRead(system.file("extdata", "test1.myCpG.txt",
                                     package = "methylKit"), sample.id="test",
                       assembly="hg38", treatment=1, context="CpG", mincov = 0))

## ---- message = FALSE, output.lines=10----------------------------------------
data("probe_ids")
probe_ids

## ---- output.lines=10---------------------------------------------------------
mapped_WGBS_data <- BSmeth2Probe(probe_id_locations = probe_ids, 
                                 WGBS_data = WGBS_GRanges,
                                 multipleMapping = TRUE,
                                 cutoff = 100)
head(mapped_WGBS_data)

## -----------------------------------------------------------------------------
deconvolution <- deconvolute(reference = extended_matrix, 
                             bulk = mapped_WGBS_data)
deconvolution[[1]]

## -----------------------------------------------------------------------------
library(granulator)
#To load the data from granulator package
load_ABIS()

#Read the bulk RNA-seq data
bulk_RNA <- bulkRNAseq_ABIS[1:100,] %>% 
  as.data.frame() %>% 
  mutate(IDs = rownames(bulkRNAseq_ABIS[1:100,])) %>%
  select("IDs", everything())

head(bulk_RNA[,1:5])

## -----------------------------------------------------------------------------
#Read the reference RNAseq data
reference_RNA <- sigMatrix_ABIS_S0 %>%
  as.data.frame() %>% 
  mutate(IDs = rownames(sigMatrix_ABIS_S0))%>%
  select("IDs", everything())

head(reference_RNA[1:5])

## -----------------------------------------------------------------------------
deconv_RNA <- deconvR::deconvolute(reference = reference_RNA,
                             bulk = bulk_RNA,model = "qp")

## -----------------------------------------------------------------------------
deconv_RNA[[1]]

## -----------------------------------------------------------------------------
sessionInfo()
stopCluster(cl)

