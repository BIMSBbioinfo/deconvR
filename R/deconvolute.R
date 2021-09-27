#' A function to deconvolute of bulk samples to their origin proportions
#' using data from reference atlas (e.g. methylation signatures)
#' Results of model are returned in dataframe "results".
#' Summary of partial R-squared values of model (min, median, mean, max...) are
#' printed upon completion.
#' @importFrom foreach %dopar%
#' @param reference A dataframe containing signatures of different cell types
#' (e.g. methylation signature) used to train the model. The first column
#' should contain a unique ID (e.g. target ID) to match rows of the reference
#' to rows of the bulk. All subsequent columns are cell types. One row per unit
#' of the signature (e.g. CpG). Each cell contains the value for the cell type
#' of this unit (e.g. methylation value of the CpG). If not given, defaults to
#' a reference atlas which is included in this package.
#' (see deconvR/inst/reference_atlas_nodup.RDS).
#' This reference atlas comes from Moss et al. (2018)
#' @param vec The user may provide a vector with which partial R-squared of the
#' results will be calculated.
#' The length must match the number of rows of the reference and bulk tables
#' merged on the ID column (with NAs removed).
#' Defaults to row means of reference.
#' @param bulk A dataframe containing signatures of bulk samples used to test
#' to model.
#' Should be dataframe with first column with unique IDs (does not need to
#' exactly match list of IDs in reference,but should have significant overlap),
#' and rest of columns = samples. Should not have duplicate IDs. May use
#' simulateCellMix function to create this dataframe.
#' @param model A string indicating which model is used to deconvolute the
#' samples. Can be either "nnls" (for non-negative least squares) or "svr"
#' (support vector regression) or "qp" (quadratic programming) or "rlm" (robust
#' linear regression). If not given, defaults to "nnls".
#' @keywords deconvolution
#' @examples
#' results_nnls <- deconvolute(bulk = simulateCellMix(10)[[1]])
#' results_rlm <- deconvolute(
#'     reference = readRDS(system.file("reference_atlas_nodup.RDS",
#'         package = "deconvR"
#'     )),
#'     bulk = simulateCellMix(5)[[1]], model = "rlm"
#' )
#' @return A list, first is a dataframe which contains predicted cell-type
#' proportions of bulk methylation profiles in "bulk", second is a list of
#' partial-rsq values of results, one value per sample.
#' @references Moss, J. et al.  (2018). Comprehensive human cell-type
#' methylation atlas reveals origins of circulating cell-free DNA in health
#' and disease. Nature communications, 9(1), 1-12.
#' \url{https://doi.org/10.1038/s41467-018-07466-6}
#' @export

deconvolute <- function(reference =
        readRDS(system.file("reference_atlas_nodup.RDS",
            package = "deconvR"
        )),
    vec = NULL, bulk, model = "nnls") {
    message("DECONVOLUTION WITH ", toupper(model))

    comb <- function(x, ...) {
        lapply(
            seq_along(x),
            function(i) c(x[[i]], lapply(list(...), function(y) y[[i]]))
        )
    }

    h <- NULL # assign a number so you can reuse it
    bulk <- tidyr::drop_na(bulk)
    ## get rid of rows in both tables with na values
    reference <- tidyr::drop_na(reference)
    find_partial_rsq <- function(observed, predicted, ref, vec) {
        if (is.null(vec)) {
            vec <- rowMeans(ref)
            ## vector defaults to row means of reference
        } else {
            assertthat::assert_that(length(vec) == length(predicted),
                msg = paste(
                    "vector should be length",
                    length(predicted), "but is length",
                    length(vec)
                )
            )
        }

        rsq_partial <- rsq::rsq.partial(
            stats::lm(predicted ~ observed + vec),
            stats::lm(predicted ~ vec)
        )$partial.rsq
        # calculate partial r-squared
        return(rsq_partial)
    }

    train_model <- function(model, ref, mix) {
        if (model == "nnls") {
            model <- nnls::nnls(data.matrix(ref), mix)
            coefficients <- model$x
        } else if (model == "svr") { # support vector regression
            model <- e1071::best.tune("svm",
                train.x = mix ~ ref, kernel = "linear",
                type = "nu-regression", scale = FALSE,
                ranges = list(nu = seq(0.25, 0.5, 0.75))
            )
            coefficients <- t(model$coefs) %*% model$SV
            coefficients <- ifelse(coefficients < 0, 0, coefficients)
        } else if (model == "qp") { # quadratic programming
            beq <- c(1)
            bvec <- c(beq, rep(0, ncol(reference[, -1])))
            Aeq <- matrix(rep(1, ncol(reference[, -1])), nrow = 1)
            Amat <- rbind(Aeq, diag(ncol(reference[, -1])))
            meq <- 1
            Dmat <- t(ref) %*% ref
            dvec <- t(ref) %*% mix

            coefficients <- quadprog::solve.QP(
                Dmat = Dmat, dvec = dvec,
                Amat = t(Amat), bvec = bvec,
                meq = meq
            )$solution

            coefficients <- ifelse(coefficients < 0, 0, coefficients)
        } else if (model == "rlm") { # robust linear regression
            model <- suppressWarnings(MASS::rlm(ref, mix, maxit = 100))
            # rlm via re-weighted least squares with maximum 100 iterations
            coefficients <- model$coefficients
            coefficients <- ifelse(coefficients < 0, 0, coefficients)
        }
    }
    # non negative least squares
    oper <- foreach::foreach(
        h = seq(2, ncol(bulk)), .inorder = TRUE,
        .combine = "comb", .multicombine = TRUE,
        .init = list(c(), list())
    ) %dopar% {
        # it's implemented this way to work with doParallel
        thedata <- tidyr::drop_na(merge(dplyr::select(bulk, 1, h),
            reference,
            by = "IDs"
        ))[, -1]
        # merge each sample to reference table,
        mix <- as.matrix(thedata[, 1])
        # first column is the bulk mixed sample
        ref <- as.matrix(thedata[, -1])
        # the rest of columns come from reference
        coefficients <- train_model(model, ref, mix)
        if (model == "svr") {
            sumOfCof <- sum(coefficients)
            coefficients <- as.numeric(coefficients / sumOfCof)
        } else {
            sumOfCof <- sum(coefficients)
            coefficients <- coefficients / sumOfCof
        }
        return(list(
            find_partial_rsq(mix, ref %*% coefficients, ref, vec),
            coefficients
        ))
    }

    rsq_partial <- oper[[1]]
    # results table will have coefficient predictions of each sample
    get_res <- function(i) {
        res <- unlist(oper[[2]][[i]])
        return(res)
    }
    expect <- NCOL(reference) -1
    results <- data.frame(t(
        vapply(seq_along(oper[[2]]), get_res, FUN.VALUE = numeric(expect))
    ))
    rownames(results) <- colnames(bulk)[-1]
    colnames(results) <- colnames(reference)[-1]
    message("SUMMARY OF PARTIAL R-SQUARED VALUES FOR ", toupper(model), ": ")
    print(summary(unlist(rsq_partial)))

    return(list(results, rsq_partial))
}
