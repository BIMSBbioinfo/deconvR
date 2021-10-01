#' @title A function to extract partial R-squares from the model results based
#' on users desicions
#' @importFrom assertthat assert_that
#' @importFrom rsq rsq.partial
#' @importFrom stats lm
#' @param vec The user may provide a vector with which partial R-squared of the
#' results will be calculated.
#' The length must match the number of rows of the reference and bulk tables
#' merged on the ID column (with NAs removed).
#' Defaults to row means of reference.
#' @param observed the observed value of the dependent variable
#' @param predicted the predicted value of the dependent variable
#' @param ref reference atlas provided by the user in matrix form
#' @param mix bulk data provided by the user in matrix form
#' @return a partial Rsuare value obtained from the result of the model
#' @keywords internal
#' @noRd
findPartialRsquare <- function(observed, predicted, ref, vec) {
    ## vector defaults to row means of reference
    if (is.null(vec)) {
        vec <- rowMeans(ref)
    } else {
        assert_that(length(vec) == length(predicted),
            msg = paste(
                "vector should be length",
                length(predicted), "but is length",
                length(vec)
            )
        )
    }
    # calculate partial r-squared
    rsq_partial <- rsq.partial(
        lm(predicted ~ observed + vec),
        lm(predicted ~ vec)
    )$partial.rsq
}
