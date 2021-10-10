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

## ---- message = FALSE, output.lines=10----------------------------------------
data("IlluminaMethEpicB5ProbeIDs")
head(IlluminaMethEpicB5ProbeIDs)

## ---- message = FALSE, output.lines=10----------------------------------------
library(deconvR) 

data("HumanCellTypeMethAtlas")
head(HumanCellTypeMethAtlas)

## ---- message = FALSE, output.lines=10----------------------------------------
samples <- simulateCellMix(3,reference = HumanCellTypeMethAtlas)[[1]]
head(samples)

## ---- message = FALSE, output.lines=10----------------------------------------
sampleMeta <- data.table("Experiment_accession" = colnames(samples)[-1],
                         "Biosample_term_name" = "new cell type")
head(sampleMeta)

## ---- output.lines=10---------------------------------------------------------
extended_matrix <- findSignatures(samples = samples, 
                                 sampleMeta = sampleMeta, 
                                 atlas = HumanCellTypeMethAtlas)
head(extended_matrix)

## ---- message = FALSE, output.lines=10----------------------------------------
WGBS_GRanges<- readRDS(system.file("extdata", "WGBS_GRanges.RDS",
                                     package = "deconvR"))
WGBS_GRanges

## ---- message = FALSE, output.lines=10----------------------------------------
WGBS_data <- readRDS(system.file("extdata", "WGBS_methylkit.RDS",
                                     package = "deconvR"))
head(WGBS_data)

## ---- message = FALSE, output.lines=10----------------------------------------
data("IlluminaMethEpicB5ProbeIDs")
head(IlluminaMethEpicB5ProbeIDs)

## ---- output.lines=10---------------------------------------------------------
mapped_WGBS_data <- BSmeth2Probe(probe_id_locations = IlluminaMethEpicB5ProbeIDs, 
                                 WGBS_data = WGBS_GRanges,
                                 multipleMapping = TRUE,
                                 cutoff = 10)
head(mapped_WGBS_data)

## -----------------------------------------------------------------------------
deconvolution <- deconvolute(reference = HumanCellTypeMethAtlas, 
                             bulk = mapped_WGBS_data)
deconvolution[[1]]

## -----------------------------------------------------------------------------
library(granulator)
#To load the data from granulator package
load_ABIS()

#Read the bulk RNA-seq data
bulk_RNA <- bulkRNAseq_ABIS[1:50,] %>% 
  as.data.frame() %>% 
  mutate(IDs = rownames(bulkRNAseq_ABIS[1:50,])) %>%
  select("IDs", everything())

head(bulk_RNA[,1:5])

## -----------------------------------------------------------------------------
#Read the reference RNAseq data
reference_RNA <- sigMatrix_ABIS_S0 %>%
  as.data.frame() %>% 
  mutate(IDs = rownames(sigMatrix_ABIS_S0))%>%
  select("IDs", everything())

head(reference_RNA[1:5])

## ---- message = FALSE, output.lines=10----------------------------------------
samples <- simulateCellMix(3,reference = reference_RNA)[[1]]
head(samples)

## ---- message = FALSE, output.lines=10----------------------------------------
sampleMeta <- data.table("Experiment_accession" = colnames(samples)[-1],
                         "Biosample_term_name" = "new cell type")
head(sampleMeta)

## ---- output.lines=10---------------------------------------------------------
extended_matrix <- findSignatures(samples = samples, 
                                  sampleMeta = sampleMeta, 
                                  atlas = reference_RNA)
head(extended_matrix)

## -----------------------------------------------------------------------------
deconv_RNA <- deconvR::deconvolute(reference = reference_RNA,
                             bulk = bulk_RNA,model = "qp")

## -----------------------------------------------------------------------------
head(deconv_RNA[[1]])

## -----------------------------------------------------------------------------
sessionInfo()
stopCluster(cl)

