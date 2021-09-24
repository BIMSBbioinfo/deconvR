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
library(ggplot2)
library(reshape2)
library(doParallel)
library(dplyr)

cl <- parallel::makeCluster(2)
doParallel::registerDoParallel(cl)
listboth100 = simulateCellMix(100)
propstable100 = listboth100[[2]]
bulktable100 = listboth100[[1]]
bulktable10 = simulateCellMix(10)[[1]]
bulktable1 = simulateCellMix(1)[[1]]

## ----eval=FALSE---------------------------------------------------------------
#  devtools::install_github("BIMSBbioinfo/deconvR")

## ----eval=FALSE---------------------------------------------------------------
#  simulateCellMix(
#    numberOfSamples,
#    mixingVector=NULL,
#    reference = readRDS(system.file("reference_atlas_nodup.RDS",
#                                     package = "deconvR")))

## ----eval=FALSE---------------------------------------------------------------
#  deconvolute(
#    reference = readRDS(system.file("reference_atlas_nodup.RDS",
#                                    package = "deconvR")),vec = NULL, bulk,
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
atlas <- readRDS(system.file("reference_atlas_nodup.RDS", package = "deconvR"))
head(atlas)

## ---- message = FALSE, output.lines=10----------------------------------------
WGBS_data <- readRDS(system.file("WGBS_GRanges.RDS", package = "deconvR"))
WGBS_data

## ---- message = FALSE, output.lines=10----------------------------------------
probe_id_locations <- readRDS(system.file("illumina_probes_hg38_GRanges.RDS", package = "deconvR"))
probe_id_locations

## ---- message = FALSE, output.lines=10----------------------------------------
atlas <- readRDS(system.file("reference_atlas_nodup.RDS", package = "deconvR"))
head(atlas)

## ---- message = FALSE, output.lines=10----------------------------------------
samples <- simulateCellMix(3)[[1]]
head(samples)

## ---- message = FALSE, output.lines=10----------------------------------------
sampleMeta <- data.table("Experiment_accession" = colnames(samples)[-1],
                         "Biosample_term_name" = "new cell type")
head(sampleMeta)

## ---- output.lines=10---------------------------------------------------------
extended_matrix <- findSignatures(samples = samples, 
                                 sampleMeta = sampleMeta, 
                                 atlas = atlas)
head(extended_matrix)

## ---- message = FALSE, output.lines=10----------------------------------------
WGBS_data <- readRDS(system.file("WGBS_GRanges.RDS", package = "deconvR"))
WGBS_data

## ---- message = FALSE, output.lines=10----------------------------------------
head(methylKit::methRead(system.file("extdata", "test1.myCpG.txt",
                                     package = "methylKit"), sample.id="test",
                       assembly="hg38", treatment=1, context="CpG", mincov = 0))

## ---- message = FALSE, output.lines=10----------------------------------------
probe_id_locations <- readRDS(system.file("illumina_probes_hg38_GRanges.RDS",
                                          package = "deconvR"))
probe_id_locations

## ---- output.lines=10---------------------------------------------------------
mapped_WGBS_data <- BSmeth2Probe(probe_id_locations = probe_id_locations, 
                                 WGBS_data = WGBS_data,
                                 multipleMapping = TRUE,
                                 cutoff = 100)
head(mapped_WGBS_data)

## ---- output.lines=10---------------------------------------------------------
deconvolution <- deconvolute(reference = extended_matrix, 
                             bulk = mapped_WGBS_data)
deconvolution[[1]]

## -----------------------------------------------------------------------------
nnls_result_100 <- deconvolute(bulk=bulktable100, model="nnls")[[1]]

## -----------------------------------------------------------------------------
svr_result_100 <- deconvolute(bulk=bulktable100, model="svr")[[1]]

## -----------------------------------------------------------------------------
qp_result_100 <- deconvolute(bulk=bulktable100, model="qp")[[1]]

## -----------------------------------------------------------------------------
rlm_result_100 <- deconvolute(bulk=bulktable100, model="rlm")[[1]]

## -----------------------------------------------------------------------------
nnls_residuals_100 <- nnls_result_100 - propstable100
nnls_rmse <- sqrt(rowMeans(nnls_residuals_100^2))
nnls_rmse <- data.frame(nnls_rmse)
summary(nnls_rmse)

## -----------------------------------------------------------------------------
boxplot(nnls_rmse)

## -----------------------------------------------------------------------------
svr_residuals_100 <- svr_result_100 - propstable100
svr_rmse <- sqrt(rowMeans(svr_residuals_100^2))
svr_rmse <- data.frame(svr_rmse)

summary(svr_rmse)

## -----------------------------------------------------------------------------
boxplot(svr_rmse)

## -----------------------------------------------------------------------------
qp_residuals_100 <- qp_result_100 - propstable100
qp_rmse <- sqrt(rowMeans(qp_residuals_100^2))
qp_rmse <- data.frame(qp_rmse)

summary(qp_rmse)

## -----------------------------------------------------------------------------
boxplot(qp_rmse)

## -----------------------------------------------------------------------------
rlm_residuals_100 <- rlm_result_100 - propstable100
rlm_rmse <- sqrt(rowMeans(rlm_residuals_100^2))
rlm_rmse <- data.frame(rlm_rmse)

summary(rlm_rmse)

## -----------------------------------------------------------------------------
boxplot(rlm_rmse)

## ---- fig.show='hide'---------------------------------------------------------
nnls_time1 <- system.time(deconvolute(bulk=bulktable1, model = "nnls"))[3]
nnls_time10 <- system.time(deconvolute(bulk=bulktable10, model = "nnls"))[3]
nnls_time100 <- system.time(deconvolute(bulk=bulktable100, model = "nnls"))[3]

print(paste("Time To Deconvolute 1 Sample:", nnls_time1,"seconds (",nnls_time1/1,"seconds/sample )"))
print(paste("Time To Deconvolute 10 Samples:", nnls_time10,"seconds (",nnls_time10/10,"seconds/sample )"))
print(paste("Time To Deconvolute 100 Samples:", nnls_time100,"seconds (",nnls_time100/100,"seconds/sample )"))



## ---- fig.show='hide'---------------------------------------------------------
svr_time1 <- system.time(deconvolute(bulk=bulktable1, model = "svr"))[3]
svr_time10 <- system.time(deconvolute(bulk=bulktable10, model = "svr"))[3]
svr_time100 <- system.time(deconvolute(bulk=bulktable100, model = "svr"))[3]

print(paste("Time To Deconvolute 1 Sample:", svr_time1,"seconds (",svr_time1/1,"seconds/sample)"))
print(paste("Time To Deconvolute 10 Samples:", svr_time10,"seconds (",svr_time10/10,"seconds/sample)"))
print(paste("Time To Deconvolute 100 Samples:", svr_time100,"seconds (",svr_time100/100,"seconds/sample)"))

## ---- fig.show='hide'---------------------------------------------------------
qp_time1 <- system.time(deconvolute(bulk=bulktable1, model = "qp"))[3]
qp_time10 <- system.time(deconvolute(bulk=bulktable10, model = "qp"))[3]
qp_time100 <- system.time(deconvolute(bulk=bulktable100, model = "qp"))[3]

print(paste("Time To Deconvolute 1 Sample:", qp_time1,"seconds (",qp_time1/1,"seconds/sample)"))
print(paste("Time To Deconvolute 10 Samples:", qp_time10,"seconds (",qp_time10/10,"seconds/sample)"))
print(paste("Time To Deconvolute 100 Samples:", qp_time100,"seconds (",qp_time100/100,"seconds/sample)"))



## ---- fig.show='hide'---------------------------------------------------------
rlm_time1 <- system.time(deconvolute(bulk=bulktable1, model = "rlm"))[3]
rlm_time10 <- system.time(deconvolute(bulk=bulktable10, model = "rlm"))[3]
rlm_time100 <- system.time(deconvolute(bulk=bulktable100, model = "rlm"))[3]

print(paste("Time To Deconvolute 1 Sample:", rlm_time1,"seconds (",rlm_time1/1,"seconds/sample)"))
print(paste("Time To Deconvolute 10 Samples:", rlm_time10,"seconds (",rlm_time10/10,"seconds/sample)"))
print(paste("Time To Deconvolute 100 Samples:", rlm_time100,"seconds (",rlm_time100/100,"seconds/sample)"))



## -----------------------------------------------------------------------------
time_data <- data.frame(c(1,10,100), 
                        c(nnls_time1, nnls_time10, nnls_time100),
                        c(svr_time1, svr_time10, svr_time100),
                        c(qp_time1, qp_time10, qp_time100),
                        c(rlm_time1, rlm_time10, rlm_time100))
colnames(time_data) <- c("NumberOfSamples", "NNLS", "SVR", "QP", "RLM")

ggplot2::ggplot(time_data, ggplot2::aes(x=NumberOfSamples)) + 
  ggplot2::geom_line(ggplot2::aes(y = NNLS, color = "steelblue")) + 
  ggplot2::geom_line(ggplot2::aes(y = SVR, color="green")) +
  ggplot2::geom_line(ggplot2::aes(y = QP, color = "pink"),linetype="twodash") +
  ggplot2::geom_line(ggplot2::aes(y = RLM, color="orange"),linetype="twodash") +
  ylab("Time (seconds)")+  xlab("Number of Samples")+
  scale_color_discrete(name = "Model", labels = c("SVR", "RLM","QP", "NNLS"))+
  ggplot2::ggtitle("Runtime of Models")



## -----------------------------------------------------------------------------
samplenumslist <- c()
for (i in 1:100) {samplenumslist <- append(samplenumslist, paste("Sample",i))}
scatterdf_nnls <- data.frame(propstable100[,1],nnls_result_100[,1], 
                             colnames(propstable100)[1])
colnames(scatterdf_nnls) <- c("actual", "predicted", "celltype")
for (i in 2:25) {
    scatterdf1 <- data.frame(propstable100[,i], nnls_result_100[,i], 
                             colnames(propstable100)[i])
    colnames(scatterdf1) <- c("actual", "predicted", "celltype")
    scatterdf_nnls <- rbind(scatterdf_nnls, scatterdf1)
}
scatterdf_nnls[,4] <- samplenumslist
colnames(scatterdf_nnls)[4] <- "samplenum"
sp_nnls <- ggplot2::ggplot(scatterdf_nnls, 
                           ggplot2::aes(x=actual, y=predicted,color=celltype)) +
  ggplot2::geom_point() + 
  ggplot2::ggtitle("NNLS Overall")
sp_nnls

## -----------------------------------------------------------------------------
nnls_bestcase <- subset(nnls_rmse, 
                        nnls_rmse <= quantile(nnls_rmse, 
                                              probs = c(0, 0.05, 0.95, 1))[2])
scatter_nnls_bestcase <- subset(scatterdf_nnls, 
                                is.element(scatterdf_nnls[,4],
                                           rownames(nnls_bestcase)))
sp_nnls_bestcase <- ggplot2::ggplot(scatter_nnls_bestcase, 
                                    ggplot2::aes(x=actual, y=predicted,
                                                 color=celltype)) + ggplot2::geom_point()+ 
  ggplot2::ggtitle("NNLS Best Case")
sp_nnls_bestcase

## -----------------------------------------------------------------------------
nnls_worstcase <- subset(nnls_rmse, 
                         nnls_rmse >= quantile(nnls_rmse,
                                               probs = c(0, 0.05, 0.95, 1))[3])
scatter_nnls_worstcase <- subset(scatterdf_nnls, 
                                 is.element(scatterdf_nnls[,4],
                                            rownames(nnls_worstcase)))
sp_nnls_worstcase <- ggplot2::ggplot(scatter_nnls_worstcase,
                                     ggplot2::aes(x=actual, y=predicted,
                                                  color=celltype)) +
  ggplot2::geom_point()+ 
  ggplot2::ggtitle("NNLS Worst Case")
sp_nnls_worstcase

## -----------------------------------------------------------------------------
scatterdf_svr <- data.frame(propstable100[,1], svr_result_100[,1], 
                            colnames(propstable100)[1])
colnames(scatterdf_svr) <- c("actual", "predicted", "celltype")
for (i in 2:25) {
    scatterdf1 <- data.frame(propstable100[,i], svr_result_100[,i], 
                             colnames(propstable100)[i])
    colnames(scatterdf1) <- c("actual", "predicted", "celltype")
    scatterdf_svr <- rbind(scatterdf_svr, scatterdf1)
}
scatterdf_svr[,4] <- samplenumslist
colnames(scatterdf_svr)[4] <- "samplenum"
sp_svr <- ggplot2::ggplot(scatterdf_svr, 
                          ggplot2::aes(x=actual, y=predicted, color=celltype))+
  ggplot2::geom_point() +
  ggplot2::ggtitle("SVR Ovrerll")
sp_svr

## -----------------------------------------------------------------------------
svr_bestcase <- subset(svr_rmse, 
                       svr_rmse <= quantile(svr_rmse, 
                                            probs = c(0, 0.05, 0.95, 1))[2])
scatter_svr_bestcase <- subset(scatterdf_svr, is.element(scatterdf_svr[,4],
                                                         rownames(svr_bestcase)))
sp_svr_bestcase <- ggplot2::ggplot(scatter_svr_bestcase, 
                                   ggplot2::aes(x=actual, y=predicted, 
                                                color=celltype)) +
  ggplot2::geom_point()+
  ggplot2::ggtitle("SVR Best Case")
sp_svr_bestcase

## -----------------------------------------------------------------------------
svr_worstcase <- subset(svr_rmse, 
                        svr_rmse >= quantile(svr_rmse, 
                                             probs = c(0, 0.05, 0.95, 1))[3])
scatter_svr_worstcase <- subset(scatterdf_svr, 
                                is.element(scatterdf_svr[,4], 
                                           rownames(svr_worstcase)))
sp_svr_worstcase <- ggplot2::ggplot(scatter_svr_worstcase, 
                                    ggplot2::aes(x=actual, y=predicted,
                                                 color=celltype)) +
  ggplot2::geom_point()+ 
  ggplot2::ggtitle("SVR Worst Case")
sp_svr_worstcase

## -----------------------------------------------------------------------------
scatterdf_qp <- data.frame(propstable100[,1], qp_result_100[,1], 
                           colnames(propstable100)[1])
colnames(scatterdf_qp) <- c("actual", "predicted", "celltype")
for (i in 2:25) {
    scatterdf1 <- data.frame(propstable100[,i], qp_result_100[,i], 
                             colnames(propstable100)[i])
    colnames(scatterdf1) <- c("actual", "predicted", "celltype")
    scatterdf_qp <- rbind(scatterdf_qp, scatterdf1)
}
scatterdf_qp[,4] <- samplenumslist
colnames(scatterdf_qp)[4] <- "samplenum"
sp_qp <- ggplot2::ggplot(scatterdf_qp, 
                         ggplot2::aes(x=actual, y=predicted, color=celltype)) +
  ggplot2::geom_point() +
  ggplot2::ggtitle("QP Overall")
sp_qp

## -----------------------------------------------------------------------------

qp_bestcase <- subset(qp_rmse, 
                      qp_rmse <= quantile(qp_rmse, 
                                          probs = c(0, 0.05, 0.95, 1))[2])
scatter_qp_bestcase <- subset(scatterdf_qp, is.element(scatterdf_qp[,4],
                                                       rownames(qp_bestcase)))
sp_qp_bestcase <- ggplot2::ggplot(scatter_qp_bestcase, 
                                  ggplot2::aes(x=actual, y=predicted, 
                                               color=celltype)) +
  ggplot2::geom_point()+
  ggplot2::ggtitle("QP Best Case")
sp_qp_bestcase

## -----------------------------------------------------------------------------
qp_worstcase <- subset(qp_rmse, 
                       qp_rmse >= quantile(qp_rmse, 
                                           probs = c(0, 0.05, 0.95, 1))[3])
scatter_qp_worstcase <- subset(scatterdf_qp, is.element(scatterdf_qp[,4],
                                                        rownames(qp_worstcase)))
sp_qp_worstcase <- ggplot2::ggplot(scatter_qp_worstcase, 
                                   ggplot2::aes(x=actual, y=predicted,
                                                color=celltype)) +
  ggplot2::geom_point()+
  ggplot2::ggtitle("QP Worst Case")
sp_qp_worstcase

## -----------------------------------------------------------------------------
scatterdf_rlm <- data.frame(propstable100[,1], rlm_result_100[,1],
                            colnames(propstable100)[1])
colnames(scatterdf_rlm) <- c("actual", "predicted", "celltype")
for (i in 2:25) {
    scatterdf1 <- data.frame(propstable100[,i], rlm_result_100[,i],
                             colnames(propstable100)[i])
    colnames(scatterdf1) <- c("actual", "predicted", "celltype")
    scatterdf_rlm <- rbind(scatterdf_rlm, scatterdf1)
}
scatterdf_rlm[,4] <- samplenumslist
colnames(scatterdf_rlm)[4] <- "samplenum"
sp_rlm <- ggplot2::ggplot(scatterdf_rlm, 
                          ggplot2::aes(x=actual, y=predicted, color=celltype)) +
  ggplot2::geom_point() +
  ggplot2::ggtitle("RLM Overall")
sp_rlm

## -----------------------------------------------------------------------------

rlm_bestcase <- subset(rlm_rmse, 
                       rlm_rmse <= quantile(rlm_rmse, 
                                            probs = c(0, 0.05, 0.95, 1))[2])
scatter_rlm_bestcase <- subset(scatterdf_rlm, 
                               is.element(scatterdf_rlm[,4],
                                          rownames(rlm_bestcase)))
sp_rlm_bestcase <- ggplot2::ggplot(scatter_rlm_bestcase, 
                                   ggplot2::aes(x=actual, y=predicted,
                                                color=celltype)) +
  ggplot2::geom_point()+ggplot2::ggtitle("RLM Best Case")
sp_rlm_bestcase

## -----------------------------------------------------------------------------
rlm_worstcase <- subset(rlm_rmse, 
                        rlm_rmse >= quantile(rlm_rmse, 
                                             probs = c(0, 0.05, 0.95, 1))[3])
scatter_rlm_worstcase <- subset(scatterdf_rlm, 
                                is.element(scatterdf_rlm[,4],
                                           rownames(rlm_worstcase)))
sp_rlm_worstcase <- ggplot2::ggplot(scatter_rlm_worstcase, 
                                    ggplot2::aes(x=actual, y=predicted,
                                                 color=celltype)) +
  ggplot2::geom_point()+ggplot2::ggtitle("RLM Worst Case")
sp_rlm_worstcase

## -----------------------------------------------------------------------------
nnls_absolute_residuals <- abs(nnls_residuals_100)
melt_nnls <- reshape2::melt(nnls_absolute_residuals)
plt_nnls <- ggplot2::ggplot(data = melt_nnls, 
                            ggplot2::aes(x = variable, y = value))
plt_nnls + ggplot2::geom_boxplot() + ggplot2::theme_minimal() + 
  ggplot2::labs(x = "Cell Type", y = "Absolute Error") +
  ggplot2::theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ggplot2::scale_y_continuous(expand=expand_scale(mult=c(0,0.1))) +
  ggplot2::ggtitle("NNLS Absolute Error by Cell Type")

## -----------------------------------------------------------------------------
svr_absolute_residuals <- abs(svr_residuals_100)
melt_svr <- reshape2::melt(svr_absolute_residuals)
plt_svr <- ggplot2::ggplot(data = melt_svr, 
                           ggplot2::aes(x = variable,y = value))
plt_svr + ggplot2::geom_boxplot() + ggplot2::theme_minimal() + 
  ggplot2::labs(x = "Cell Type", y = "Absolute Error") +
  ggplot2::theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ggplot2::scale_y_continuous(expand=expand_scale(mult=c(0,0.1))) +
  ggplot2::ggtitle("SVR Absolute Error by Cell Type")

## -----------------------------------------------------------------------------
qp_absolute_residuals <- abs(qp_residuals_100)
melt_qp <- reshape2::melt(qp_absolute_residuals)
plt_qp <- ggplot2::ggplot(data = melt_qp , 
                          ggplot2::aes(x = variable,y = value))
plt_qp + ggplot2::geom_boxplot() + ggplot2::theme_minimal() + 
  ggplot2::labs(x = "Cell Type", y = "Absolute Error") +
  ggplot2::theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ggplot2::scale_y_continuous(expand=expand_scale(mult=c(0,0.1))) +
  ggplot2::ggtitle("QP Absolute Error by Cell Type")

## -----------------------------------------------------------------------------
rlm_absolute_residuals <- abs(rlm_residuals_100)
melt_rlm <- reshape2::melt(rlm_absolute_residuals)
plt_rlm <- ggplot2::ggplot(data = melt_rlm , 
                           ggplot2::aes(x = variable, y = value))
plt_rlm + ggplot2::geom_boxplot() + ggplot2::theme_minimal() + 
  ggplot2::labs(x = "Cell Type", y = "Absolute Error") +
  ggplot2::theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ggplot2::scale_y_continuous(expand=expand_scale(mult=c(0,0.1))) + 
  ggplot2::ggtitle("RLM Absolute Error by Cell Type")

## -----------------------------------------------------------------------------
sessionInfo()
stopCluster(cl)

