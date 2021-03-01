#' A function to deconvolute of cell-free DNA methylation signals to their cell type of origin.
#' Results of model are returned in dataframe "results".
#' Summary of RMSE values of model (min, median, mean, max...) are printed upon completion.
#' @importFrom foreach %dopar%
#' @param reference A dataframe containing CpG signatures of different cell types used to train the model. The first column
#' should contain a unique ID (e.g. target ID) to match rows of the reference to rows of the bulk. All subsequent columns are cell types.
#'  Rows are CpGs. Each cell contains the methylation value for the cell type at the CpG location. If not given, defaults to a reference atlas
#'  which is included in this package (see deconvR/inst/reference_atlas_nodup.RDS). This atlas contains 25 cell types, 6105 CpGs (and their target IDs),
#'  and therefore 152,625 methylation values.
#' @param bulk A dataframe containing CpG signatures of different bulk methylation profiles used to test to model.
#' Should be dataframe with first column "CpGs" for Illumina IDs (does not need to exactly match list of IDs in reference,
#' but should have significant overlap), and rest of columns = samples. Should not have duplicate Illumina IDs. May use
#' simulateCellMix function to create this dataframe.
#' @param model A string indicating which model is used to deconvolute the samples. Can be either "nnls" (for non-negative least squares) or
#' "svr" (support vector regression) or "qp" (quadratic programming) or "rlm" (robust linear regression). If not given, defaults to "nnls".
#' @keywords deconvolution
#' @examples
#'  deconvolute(bulk=simulateCellMix(50)[[1]])
#'  deconvolute(reference=readRDS(system.file("reference_atlas_nodup.RDS", package = "deconvR")),
#'     bulk=simulateCellMix(5)[[1]],model="rlm")
#' @return A dataframe which contains predicted cell-type proportions of bulk methylation profiles in "bulk".
#' @export



#options(scipen = 999)
#reference_atlas = read.csv(file = './original_reference_atlas.csv', header = TRUE, row.names = 1)
#bulk = read.csv(file = './bulk.csv', header = TRUE, row.names = 1)
#reference_atlas_nodup <- reference_atlas[!duplicated(reference_atlas$CpGs),]
#bulk_nodup <- bulk[!duplicated(bulk$CpGs),]

#install.packages("tidyverse")
#install.packages("e1071")
#install.packages("quadprog")
#install.packages("nnls")
#install.packages("MASS")
#library(tidyverse)
#library(e1071)
#library(quadprog)
#library(nnls)
#library(MASS)

deconvolute = function(reference=readRDS(system.file("reference_atlas_nodup.RDS", package = "deconvR")), bulk, model="nnls") {
  print(paste("DECONVOLUTION WITH",toupper(model)))

  #results_RMSEs = data.frame(matrix(nrow=length(colnames(bulk))-1,ncol=1, dimnames=list(colnames(bulk)[-1], c("RMSE"))))
  results = data.frame(matrix(ncol=length(colnames(reference))-1,nrow=length(colnames(bulk))-1, dimnames=list( colnames(bulk)[-1], colnames(reference)[-1])))
  comb <- function(x, ...) {
    lapply(seq_along(x),
           function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
  }
  if (model == "nnls") {    #non negative least squares
    oper <- foreach::foreach (h = 2:ncol(bulk), .inorder = TRUE, .combine = "comb", .multicombine=TRUE, .init=list(c(), c(), list())) %dopar% { #skip first column because it's CpGs
      thedata = tidyr::drop_na(merge(dplyr::select(bulk, 1, h), reference, by="CpGs"))[,-1] #merge each sample to reference table, removing CpGs that aren't present in both (NAs)
      ref_nnls = as.matrix(thedata[,-1])
      mix_nnls = as.matrix(thedata[,1])
      nnls_model = nnls::nnls(data.matrix(ref_nnls), mix_nnls)
      nnls_coefficients = nnls_model$x
      sumOfCof = sum(nnls_coefficients)
      nnls_coefficients = nnls_coefficients / sumOfCof #normalize so coefficients add to 1

      nnls_residuals = nnls_model$residuals
      #nnls_residuals = mix_nnls - nnls_model$fitted
      RSS = sum(nnls_residuals^2)
      variance = mix_nnls - mean(mix_nnls)
      TSS = sum(variance^2)
      nnls_rsq = 1 - (RSS / TSS)               #r-squared
      nnls_sumofsquares = nnls_model$deviance  #residual sum of squares


      return(list(nnls_rsq, nnls_sumofsquares, nnls_coefficients))


    }

    rsq = oper[[1]]
    rss = oper[[2]]
    for (i in 1:length(oper[[3]])){
      results[i,] = unlist(oper[[3]][[i]]) #results table will have coefficient predictions of each sample
    }


  }
  else if (model == "svr") {  #support vector regression
    oper <- foreach::foreach(h = 2:ncol(bulk), .inorder = TRUE, .combine = "comb", .multicombine=TRUE, .init=list(c(), c(), list())) %dopar% {  #skip first column because it's CpGs
      thedata = tidyr::drop_na(merge(dplyr::select(bulk, 1, h), reference, by="CpGs"))[,-1]
      ref_svr = as.matrix(thedata[,-1])
      mix_svr = as.matrix(thedata[,1])


      svr_model = e1071::best.tune("svm", train.x = mix_svr ~ ref_svr, kernel="linear", type="nu-regression", scale=FALSE, ranges = list(nu = seq(0.25,0.5,0.75))) #find the best nu values for the models

      svr_coefficients = t(svr_model$coefs) %*% svr_model$SV


      svr_coefficients= ifelse(svr_coefficients < 0, 0, svr_coefficients)

      sumOfCof = sum(svr_coefficients)
      svr_coefficients = as.numeric(svr_coefficients / sumOfCof)  #normalize so coefficients add to 1

      svr_fitted = ref_svr %*% svr_coefficients
      svr_residuals = mix_svr - svr_fitted
      RSS = sum(svr_residuals^2)
      variance = mix_svr - mean(mix_svr)
      TSS = sum(variance^2)
      svr_rsq = 1 - (RSS / TSS)               #r-squared
      svr_sumofsquares = sum(svr_residuals^2) #residual sum of squares

     return(list(svr_rsq, svr_sumofsquares, svr_coefficients))

    }
    rsq = oper[[1]]
    rss = oper[[2]]
    for (i in 1:length(oper[[3]])){
      results[i,] = unlist(oper[[3]][[i]]) #results table will have coefficient predictions of each sample
    }
  }
  else if (model == "qp") { #quadratic programming
    beq = c(1)
    bvec=c(beq, rep(0, ncol(reference[,-1])))
    Aeq = matrix(rep(1, ncol(reference[,-1])), nrow = 1)
    Amat = rbind(Aeq, diag(ncol(reference[,-1])))
    meq = 1

    oper <- foreach::foreach (h = 2:ncol(bulk), .inorder = TRUE, .combine = "comb", .multicombine=TRUE, .init=list(c(), c(), list())) %dopar% {  #skip first column because it's CpGs
      thedata = tidyr::drop_na(merge(dplyr::select(bulk, 1, h), reference, by="CpGs"))[,-1] #merge each sample to reference table, removing CpGs that aren't present in both (NAs)
      ref_qp = as.matrix(thedata[,-1])
      mix_qp = as.matrix(thedata[,1])

      Dmat = t(ref_qp) %*% ref_qp #i got these values using this tutorial https://github.com/dariober/quadratic-programming-deconvolution/blob/master/quadratic_programming_tutorial.pdf
      dvec = t(ref_qp) %*% mix_qp

      qp_coefficients = quadprog::solve.QP(Dmat = Dmat, dvec = dvec, Amat = t(Amat), bvec = bvec, meq = meq)$solution

      qp_coefficients = ifelse(qp_coefficients < 0, 0, qp_coefficients)

      sumOfCof = sum(qp_coefficients)
      qp_coefficients = qp_coefficients / sumOfCof  #normalize so coefficients add to 1

      qp_fitted = ref_qp %*% qp_coefficients
      qp_residuals = mix_qp - qp_fitted
      RSS = sum(qp_residuals^2)
      variance = mix_qp - mean(mix_qp)
      TSS = sum(variance^2)
      qp_rsq = 1 - (RSS / TSS)               #r-squared
      qp_sumofsquares = sum(qp_residuals^2)  #residual sum of squares

      return(list(qp_rsq, qp_sumofsquares, qp_coefficients))

    }

    rsq = oper[[1]]
    rss = oper[[2]]
    for (i in 1:length(oper[[3]])){
      results[i,] = unlist(oper[[3]][[i]]) #results table will have coefficient predictions of each sample
    }
  }
  else if (model == "rlm") {  #robust linear regression
    oper <- foreach::foreach (h = 2:ncol(bulk), .inorder = TRUE, .combine = "comb", .multicombine=TRUE, .init=list(c(), c(), list())) %dopar% {  #skip first column because it's CpGs
      thedata = tidyr::drop_na(merge(dplyr::select(bulk, 1, h), reference, by="CpGs"))[,-1] #merge each sample to reference table, removing CpGs that aren't present in both (NAs)
      ref_rlm = as.matrix(thedata[,-1])
      mix_rlm = as.matrix(thedata[,1])

      rlm_model = suppressWarnings(MASS::rlm(ref_rlm, mix_rlm, maxit = 100))
      rlm_coefficients = rlm_model$coefficients

      rlm_coefficients = ifelse(rlm_coefficients < 0, 0, rlm_coefficients)

      sumOfCof = sum(rlm_coefficients)

      rlm_coefficients = rlm_coefficients / sumOfCof  #normalize so coefficients add to 1

      rlm_residuals = rlm_model$residuals

      RSS = sum(rlm_residuals^2)
      variance = mix_rlm - mean(mix_rlm)
      TSS = sum(variance^2)
      rlm_rsq = 1 - (RSS / TSS)               #r-squared
      rlm_sumofsquares = sum(rlm_residuals^2) #residual sum of squares

      return(list(rlm_rsq, rlm_sumofsquares, rlm_coefficients))

    }
    rsq = oper[[1]]
    rss = oper[[2]]
    for (i in 1:length(oper[[3]])){
      results[i,] = unlist(oper[[3]][[i]]) #results table will have coefficient predictions of each sample
    }

  }



  else {
    stop("Model should be either \"nnls\" or  \"svr\" or  \"qp\" or \"rlm\"")
  }
  print(paste("SUMMARY OF R-SQUARED VALUES FOR ",toupper(model),": "))
  print(summary(unlist(rsq)))
  print(paste("SUMMARY OF RESIDUAL-SUM-OF-SQUARES VALUES FOR",toupper(model),":"))
  print(summary(unlist(rss)))
  return(results)
}
