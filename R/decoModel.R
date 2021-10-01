#' @title A function to train user selected model to deconvolute samples
#' @param model A string indicating which model is used to deconvolute the
#' samples. Can be either "nnls" (for non-negative least squares) or "svr"
#' (support vector regression) or "qp" (quadratic programming) or "rlm" (robust
#' linear regression).
#' @param ref reference atlas provided by user in matrix form
#' @param mix bulk data provided by user in matrix form
#' @importFrom nnls nnls
#' @importFrom MASS rlm
#' @importFrom e1071 best.tune
#' @importFrom quadprog solve.QP
#' @keywords internal
#' @return the coefficient of training results
#' @noRd
decoModel <- function(model, ref, mix, colnum) {
    if (model == "nnls") { # non negative least squares
        coefficients <- (nnls(data.matrix(ref), mix))$x
    } else if (model == "svr") { # support vector regression
        model <- best.tune("svm",
            train.x = mix ~ ref, kernel = "linear",
            type = "nu-regression", scale = FALSE,
            ranges = list(nu = seq(0.25, 0.5, 0.75))
        )
        coefficients <- t(model$coefs) %*% model$SV
    } else if (model == "qp") { # quadratic programming
        Amat <- rbind(
            matrix(rep(1, colnum), nrow = 1),
            diag(colnum)
        )
        coefficients <- solve.QP(
            Dmat = (t(ref) %*% ref), dvec = (t(ref) %*% mix),
            Amat = t(Amat), bvec = c(c(1), rep(0, colnum)),
            meq = 1
        )$solution
    } else if (model == "rlm") { # robust linear regression
        coefficients <- (rlm(ref, mix, maxit = 100))$coefficients
    }
    coefficients <- ifelse(coefficients < 0, 0, coefficients)
}
