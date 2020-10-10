#' A function to deconvolute of cell-free DNA methylation signals to their cell type of origin.
#' Results of model are returned in dataframe "results".
#' Summary of RMSE values of model (min, median, mean, max...) are printed upon completion.
#' @param reference A dataframe containing CpG signatures of different cell types used to train the model. Should be dataframe with first column "CpGs" for Illumina IDs, and rest of columns = cell types. Should not have duplicate Illumina IDs. If not given, defaults to ref_atlas, a reference atlas which is included in package.
#' @param bulk A dataframe containing CpG signatures of different bulk methylation profiles used to test to model. Should be dataframe with first column "CpGs" for Illumina IDs (does not need to exactly match list of IDs in reference, but should have significant overlap), and rest of columns = cell types. Should not have duplicate Illumina IDs. May use makeBigTable function to create this dataframe.
#' @param model A string indicating which model is used to deconvolute the samples. Can be either "nnls" (for non-negative least squares) or "svr" (support vector regression) or "qp" (quadratic programming) or "rlm" (robust linear regression). If not given, defaults to "nnls".
#' @keywords deconvolution
#' @examples
#'  deconvolute(bulk=makeBigTable(50))
#'  deconvolute(reference=read.csv(system.file("reference_atlas_nodup.csv", package = "deconvR")), bulk=makeBigTable(5), model="rlm" )
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

deconvolute = function(reference=read.csv(system.file("reference_atlas_nodup.csv", package = "deconvR")), bulk, model="nnls") {
  print(paste("DECONVOLUTION WITH",toupper(model)))

  results_RMSEs = data.frame(matrix(nrow=length(colnames(bulk))-1,ncol=1, dimnames=list(colnames(bulk)[-1], c("RMSE"))))
  results = data.frame(matrix(ncol=length(colnames(reference))-1,nrow=length(colnames(bulk))-1, dimnames=list( colnames(bulk)[-1], colnames(reference)[-1])))

  if (model == "nnls") {    #non negative least squares
    for (h in 2:ncol(bulk)) { #skip first column because it's CpGs
      thedata = drop_na(merge(dplyr::select(bulk, 1, h), reference, by="CpGs"))[,-1] #merge each sample to reference table, removing CpGs that aren't present in both (NAs)
      ref_nnls = as.matrix(thedata[,-1])
      mix_nnls = as.matrix(thedata[,1])
      nnls_model = nnls(data.matrix(ref_nnls), mix_nnls)
      nnls_coefficients = nnls_model$x
      sumOfCof = sum(nnls_coefficients)
      for (i in 1:length(nnls_coefficients)) {
        nnls_coefficients[i] = nnls_coefficients[i] / sumOfCof #normalize so coefficients add to 1
        results[h-1,i] = nnls_coefficients[i] #results table will have coefficient predictions of each sample
      }

      nnls_predictions = ref_nnls %*% nnls_coefficients
      nnls_residuals = mix_nnls - nnls_predictions
      nnls_RMSE = sqrt(mean(nnls_residuals^2))
      results_RMSEs[h-1,1] = nnls_RMSE #each sample will have an RMSE value for predicted values from given coefficients by model vs actual values
    }
  }
  else if (model == "svr") {  #support vector regression
    for (h in 2:ncol(bulk)) {  #skip first column because it's CpGs
      thedata = drop_na(merge(dplyr::select(bulk, 1, h), reference, by="CpGs"))[,-1]
      ref_svr = as.matrix(thedata[,-1])
      mix_svr = as.matrix(thedata[,1])

      #svr_model = svm(thedata[,1] ~ ., thedata, kernel="linear", type="nu-regression")
      tuneResult = tune(svm, mix_svr ~ ref_svr, kernel="linear", scale=FALSE, type="nu-regression", ranges = list(nu = seq(0.25,0.5,0.75))) #find the best nu values for the models

      svr_model = tuneResult$best.model #make predictions with the previously found best model
      svr_coefficients = t(svr_model$coefs) %*% svr_model$SV

      for (i in 1:length(svr_coefficients)) {
        if (svr_coefficients[i] < 0) {
          svr_coefficients[i] = 0 #normalize so all coefficients greater or equal to zero
        }
      }
      sumOfCof = sum(svr_coefficients)
      for (i in 1:length(svr_coefficients)) {
        svr_coefficients[i] = svr_coefficients[i] / sumOfCof  #normalize so coefficients add to 1
        results[h-1,i] = svr_coefficients[i] #results table will have coefficient predictions of each sample
      }

      svr_predictions = ref_svr %*% t(svr_coefficients)
      svr_residuals = mix_svr - svr_predictions
      svr_RMSE = sqrt(mean(svr_residuals^2))
      results_RMSEs[h-1,1] = svr_RMSE  #each sample will have an RMSE value for predicted values from given coefficients by model vs actual values
    }
  }
  else if (model == "qp") { #quadratic programming
    for (h in 2:ncol(bulk)) {  #skip first column because it's CpGs
      thedata = drop_na(merge(dplyr::select(bulk, 1, h), reference, by="CpGs"))[,-1] #merge each sample to reference table, removing CpGs that aren't present in both (NAs)
      ref_qp = as.matrix(thedata[,-1])
      mix_qp = as.matrix(thedata[,1])

      Dmat = t(ref_qp) %*% ref_qp #i got these values using this tutorial https://github.com/dariober/quadratic-programming-deconvolution/blob/master/quadratic_programming_tutorial.pdf
      dvec = t(ref_qp) %*% mix_qp
      Aeq = matrix(rep(1, ncol(ref_qp)), nrow = 1)
      beq = c(1)
      Amat = rbind(Aeq, diag(ncol(ref_qp)))
      bvec = c(beq, rep(0, ncol(ref_qp)))
      meq = 1

      qp_coefficients = solve.QP(Dmat = Dmat, dvec = dvec, Amat = t(Amat), bvec = bvec, meq = meq)$solution

      for (i in 1:length(qp_coefficients)) {
        if (qp_coefficients[i] < 0) {
          qp_coefficients[i] = 0  #normalize so coefficients add to 1
        }
      }
      sumOfCof = sum(qp_coefficients)
      for (i in 1:length(qp_coefficients)) {
        qp_coefficients[i] = qp_coefficients[i] / sumOfCof  #normalize so coefficients add to 1
        results[h-1,i] = qp_coefficients[i] #results table will have coefficient predictions of each sample
      }

      qp_predictions = ref_qp %*% qp_coefficients
      qp_residuals = mix_qp - qp_predictions
      qp_RMSE = sqrt(mean(qp_residuals^2))
      results_RMSEs[h-1,1] = qp_RMSE #each sample will have an RMSE value for predicted values from given coefficients by model vs actual values
    }
  }
  else if (model == "rlm") {  #robust linear regression
    for (h in 2:ncol(bulk)) {  #skip first column because it's CpGs
      thedata = drop_na(merge(dplyr::select(bulk, 1, h), reference, by="CpGs"))[,-1] #merge each sample to reference table, removing CpGs that aren't present in both (NAs)
      ref_rlm = as.matrix(thedata[,-1])
      mix_rlm = as.matrix(thedata[,1])

      rlm_model = suppressWarnings(rlm(ref_rlm, mix_rlm, maxit = 100))
      rlm_coefficients = rlm_model$coefficients
      for (i in 1:length(rlm_coefficients)) {
        if (rlm_coefficients[i] < 0) {
          rlm_coefficients[i] = 0  #normalize so coefficients add to 1
        }
      }
      sumOfCof = sum(rlm_coefficients)
      for (i in 1:length(rlm_coefficients)) {
        rlm_coefficients[i] = rlm_coefficients[i] / sumOfCof  #normalize so coefficients add to 1
        results[h-1,i] = rlm_coefficients[i] #results table will have coefficient predictions of each sample
      }

      rlm_predictions = ref_rlm %*% rlm_coefficients
      rlm_residuals = mix_rlm - rlm_predictions
      rlm_RMSE = sqrt(mean(rlm_residuals^2))
      results_RMSEs[h-1,1] = rlm_RMSE #each sample will have an RMSE value for predicted values from given coefficients by model vs actual values

    }
  }



  else {
    stop("Model should be either \"nnls\" or  \"svr\" or  \"qp\" or \"rlm\"")
  }

  print(paste("SUMMARY OF RMSE VALUES USING ",toupper(model)))
  print(summary(results_RMSEs$RMSE))
  results
}
