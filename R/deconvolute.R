#' A function to deconvolute of bulk samples to their origin cell type proportions using e.g. methylation signature
#' Results of model are returned in dataframe "results".
#' Summary of partial R-squared values of model (min, median, mean, max...) are printed upon completion.
#' @importFrom foreach %dopar%
#' @param reference A dataframe containing signatures of different cell types (e.g. methylation signature) used to train the model. The first column
#' should contain a unique ID (e.g. target ID) to match rows of the reference to rows of the bulk. All subsequent columns are cell types.
#'  One row per unit of the signature (e.g. CpG). Each cell contains the value for the cell type for this unit (e.g. methylation value of the CpG).
#'  If not given, defaults to a reference atlas which is included in this package (see deconvR/inst/reference_atlas_nodup.RDS).
#'  This reference atlas comes from Moss et al. (2018)
#' @param bulk A dataframe containing signatures of bulk samples used to test to model.
#' Should be dataframe with first column with unique IDs (does not need to exactly match list of IDs in reference,
#' but should have significant overlap), and rest of columns = samples. Should not have duplicate IDs. May use
#' simulateCellMix function to create this dataframe.
#' @param model A string indicating which model is used to deconvolute the samples. Can be either "nnls" (for non-negative least squares) or
#' "svr" (support vector regression) or "qp" (quadratic programming) or "rlm" (robust linear regression). If not given, defaults to "nnls".
#' @keywords deconvolution
#' @examples
#'  deconvolute(bulk=simulateCellMix(50)[[1]])
#'  deconvolute(reference=readRDS(system.file("reference_atlas_nodup.RDS", package = "deconvR")),
#'     bulk=simulateCellMix(5)[[1]],model="rlm")
#' @return A dataframe which contains predicted cell-type proportions of bulk methylation profiles in "bulk".
#' @references Moss, J. et al.  (2018). Comprehensive human cell-type methylation atlas reveals origins of circulating cell-free DNA in health and disease. Nature communications, 9(1), 1-12. \url{https://doi.org/10.1038/s41467-018-07466-6}
#' @export


deconvolute = function(reference=readRDS(system.file("reference_atlas_nodup.RDS", package = "deconvR")), bulk, model="nnls") {
  print(paste("DECONVOLUTION WITH",toupper(model)))

  results = data.frame(matrix(ncol=length(colnames(reference))-1,nrow=length(colnames(bulk))-1, dimnames=list( colnames(bulk)[-1], colnames(reference)[-1])))
  h = 2
  comb <- function(x, ...) {
    lapply(seq_along(x),
           function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
  }
  if (model == "nnls") {    #non negative least squares
    oper <- foreach::foreach (h = 2:ncol(bulk), .inorder = TRUE, .combine = "comb", .multicombine=TRUE, .init=list(c(), list())) %dopar% { # it's implemented this way to work with doParallel
      thedata = tidyr::drop_na(merge(dplyr::select(bulk, 1, h), reference, by="IDs"))[,-1] # merge each sample to reference table, removing IDs that aren't present in both
      mix_nnls = as.matrix(thedata[,1])   # first column is the bulk mixed sample
      ref_nnls = as.matrix(thedata[,-1])  # the rest of columns come from reference
      nnls_model = nnls::nnls(data.matrix(ref_nnls), mix_nnls)
      nnls_coefficients = nnls_model$x
      sumOfCof = sum(nnls_coefficients)
      nnls_coefficients = nnls_coefficients / sumOfCof # normalize so coefficients add to 1

      X = mix_nnls            # X = observed values
      Y = nnls_model$fitted   # Y = predicted values
      Z = rowMeans(ref_nnls)  # Z = row means of reference (after removing IDs which weren't in bulk sample)
      nnls_rsq_partial = rsq::rsq.partial(stats::lm(Y~X+Z),stats::lm(Y~Z)) # calculate partial r-squared

      return(list(nnls_rsq_partial$partial.rsq, nnls_coefficients))
    }

    rsq_partial = oper[[1]]
    for (i in 1:length(oper[[2]])){
      results[i,] = unlist(oper[[2]][[i]]) # results table will have coefficient predictions of each sample
    }
  }

  else if (model == "svr") {  # support vector regression
    oper <- foreach::foreach(h = 2:ncol(bulk), .inorder = TRUE, .combine = "comb", .multicombine=TRUE, .init=list(c(), list())) %dopar% {
      thedata = tidyr::drop_na(merge(dplyr::select(bulk, 1, h), reference, by="IDs"))[,-1] # merge each sample to reference table, removing IDs that aren't present in both
      mix_svr = as.matrix(thedata[,1])   # first column is the bulk mixed sample
      ref_svr = as.matrix(thedata[,-1])  # the rest of columns come from reference

      svr_model = e1071::best.tune("svm", train.x = mix_svr ~ ref_svr, kernel="linear", type="nu-regression", scale=FALSE, ranges = list(nu = seq(0.25,0.5,0.75))) # find the best nu values for the models

      svr_coefficients = t(svr_model$coefs) %*% svr_model$SV

      svr_coefficients= ifelse(svr_coefficients < 0, 0, svr_coefficients) # normalize so coefficients > 0

      sumOfCof = sum(svr_coefficients)
      svr_coefficients = as.numeric(svr_coefficients / sumOfCof)  # normalize so coefficients add to 1

      X = mix_svr                        # X = observed values
      Y = ref_svr %*% svr_coefficients   # Y = predicted values
      Z = rowMeans(ref_svr)              # Z = row means of reference (after removing IDs which weren't in bulk sample)
      svr_rsq_partial = rsq::rsq.partial(stats::lm(Y~X+Z),stats::lm(Y~Z)) # calculate partial r-squared

      return(list(svr_rsq_partial$partial.rsq, svr_coefficients))
    }

    rsq_partial = oper[[1]]
    for (i in 1:length(oper[[2]])){
      results[i,] = unlist(oper[[2]][[i]]) #results table will have coefficient predictions of each sample
    }
  }

  else if (model == "qp") { #quadratic programming

    # for reference about these values read solve.QP documentation
    # also helpful: https://github.com/dariober/quadratic-programming-deconvolution/blob/master/quadratic_programming_tutorial.pdf
    beq = c(1)
    bvec=c(beq, rep(0, ncol(reference[,-1])))
    Aeq = matrix(rep(1, ncol(reference[,-1])), nrow = 1)
    Amat = rbind(Aeq, diag(ncol(reference[,-1])))
    meq = 1

    oper <- foreach::foreach (h = 2:ncol(bulk), .inorder = TRUE, .combine = "comb", .multicombine=TRUE, .init=list(c(), list())) %dopar% {
      thedata = tidyr::drop_na(merge(dplyr::select(bulk, 1, h), reference, by="IDs"))[,-1] # merge each sample to reference table, removing IDs that aren't present in both
      mix_qp = as.matrix(thedata[,1])   # first column is the bulk mixed sample
      ref_qp = as.matrix(thedata[,-1])  # the rest of columns come from reference

      Dmat = t(ref_qp) %*% ref_qp
      dvec = t(ref_qp) %*% mix_qp

      qp_coefficients = quadprog::solve.QP(Dmat = Dmat, dvec = dvec, Amat = t(Amat), bvec = bvec, meq = meq)$solution

      qp_coefficients = ifelse(qp_coefficients < 0, 0, qp_coefficients) # make sure all coefficients are > 0

      sumOfCof = sum(qp_coefficients)
      qp_coefficients = qp_coefficients / sumOfCof  # make sure coefficients add to 1

      X = mix_qp                        # X = observed values
      Y = ref_qp %*% qp_coefficients    # Y = predicted values
      Z = rowMeans(ref_qp)              # Z = row means of reference (after removing IDs which weren't in bulk sample)
      qp_rsq_partial = rsq::rsq.partial(stats::lm(Y~X+Z),stats::lm(Y~Z)) # calculate partial r-squared

      return(list(qp_rsq_partial$partial.rsq, qp_coefficients))
    }

    rsq_partial = oper[[1]]
    for (i in 1:length(oper[[2]])){
      results[i,] = unlist(oper[[2]][[i]]) #results table will have coefficient predictions of each sample
    }
  }

  else if (model == "rlm") {  #robust linear regression
    oper <- foreach::foreach (h = 2:ncol(bulk), .inorder = TRUE, .combine = "comb", .multicombine=TRUE, .init=list(c(), list())) %dopar% {
      thedata = tidyr::drop_na(merge(dplyr::select(bulk, 1, h), reference, by="IDs"))[,-1] # merge each sample to reference table, removing IDs that aren't present in both
      mix_rlm = as.matrix(thedata[,1])   # first column is the bulk mixed sample
      ref_rlm = as.matrix(thedata[,-1])  # the rest of columns come from reference

      rlm_model = suppressWarnings(MASS::rlm(ref_rlm, mix_rlm, maxit = 100))  # rlm via re-weighted least squares with maximum 100 iterations
      rlm_coefficients = rlm_model$coefficients                               # (i found 100 iterations the best balence of time vs accurracy but can be changed)

      rlm_coefficients = ifelse(rlm_coefficients < 0, 0, rlm_coefficients)    # normalize coefficients so all > 0

      sumOfCof = sum(rlm_coefficients)

      rlm_coefficients = rlm_coefficients / sumOfCof   # normalize so coefficients add to 1

      X = mix_rlm                         # X = observed values
      Y = ref_rlm %*% rlm_coefficients    # Y = predicted values
      Z = rowMeans(ref_rlm)                # Z = row means of reference (after removing IDs which weren't in bulk sample)
      rlm_rsq_partial = rsq::rsq.partial(stats::lm(Y~X+Z),stats::lm(Y~Z)) # calculate partial r-squared

      return(list(rlm_rsq_partial$partial.rsq, rlm_coefficients))

    }
    rsq_partial = oper[[1]]
    for (i in 1:length(oper[[2]])){
      results[i,] = unlist(oper[[2]][[i]])  # results table will have coefficient predictions of each sample
    }

  }



  else {
    stop("Model should be either \"nnls\" or  \"svr\" or  \"qp\" or \"rlm\"")
  }
  print(paste("SUMMARY OF PARTIAL R-SQUARED VALUES FOR",toupper(model),": "))
  print(summary(unlist(rsq_partial)))

  return(results)
}
