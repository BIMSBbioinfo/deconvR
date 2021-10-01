#' @title A function to deconvolute of bulk samples to their origin proportions
#' using data from reference atlas (e.g. methylation signatures)
#' Results of model are returned in dataframe "results".
#' Summary of partial R-squared values of model (min, median, mean, max...) are
#' printed upon completion.
#' @importFrom foreach %dopar% foreach
#' @importFrom stats  na.omit
#' @importFrom dplyr select
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
#' data("reference_atlas")
#' bulk_data <- simulateCellMix(10, reference = reference_atlas)[[1]]
#' # non-least negative square regression
#' results_nnls <- deconvolute(
#'     bulk = bulk_data,
#'     reference = reference_atlas
#' )
#' # Quadric programming
#' results_qp <- deconvolute(
#'     reference = reference_atlas,
#'     bulk = bulk_data, model = "qp"
#' )
#' @return A list, first is a dataframe which contains predicted cell-type
#' proportions of bulk methylation profiles in "bulk", second is a list of
#' partial-rsq values of results, one value per sample.
#' @references Moss, J. et al.  (2018). Comprehensive human cell-type
#' methylation atlas reveals origins of circulating cell-free DNA in health
#' and disease. Nature communications, 9(1), 1-12.
#' \url{https://doi.org/10.1038/s41467-018-07466-6}
#' @export

deconvolute <- function(reference,
    vec = NULL, bulk, model = "nnls") {
    message("DECONVOLUTION WITH ", toupper(model))

    comb <- function(x, ...) {
        lapply(seq_along(x), function(i) {
            c(x[[i]], lapply(
                list(...),
                function(y) y[[i]]
            ))
        })
    }
    ## get rid of rows in both tables with na values
    clean <- lapply(list(bulk, reference), na.omit)
    colnum <- ncol(clean[[2]][, -1])
    # save internal functions & variables within foreach,use parallelization
    h <- NULL
    foreachList <- list()
    foreachList$findPartialRsquare <- findPartialRsquare
    foreachList$decoModel <- decoModel
    operation <- foreach(
        h = seq(2, ncol(clean[[1]])), .inorder = TRUE,
        .combine = "comb", .multicombine = TRUE,
        .export = c("findPartialRsquare", "decoModel"),
        .init = list(c(), list())
    ) %dopar% {
        thedata <- na.omit(merge(dplyr::select(clean[[1]], 1, h),
            clean[[2]],
            by = "IDs"
        ))[, -1]
        # merge each sample to reference table
        mix <- as.matrix(thedata[, 1])
        ref <- as.matrix(thedata[, -1])
        coefficients <- decoModel(model, ref, mix, colnum)
        if (model == "svr") {
            coefficients <- as.numeric(coefficients / sum(coefficients))
        } else {
            coefficients <- coefficients / sum(coefficients)
        }
        return(list(
            findPartialRsquare(mix,(ref %*% coefficients) , ref, vec),
            coefficients
        ))
    }

    results <- data.frame(t(
        vapply(seq_along(operation[[2]]), function(i) {
            res <- unlist(operation[[2]][[i]])
            return(res)
        }, FUN.VALUE = numeric(NCOL(clean[[2]]) - 1))
    ), row.names = colnames(clean[[1]])[-1])
    colnames(results) <- colnames(clean[[2]])[-1]
    message("SUMMARY OF PARTIAL R-SQUARED VALUES FOR ", toupper(model), ": ")
    print(summary(unlist(operation[[1]])))

    # results table will have coefficient predictions of each sample
    return(list(results, operation[[1]]))
}
