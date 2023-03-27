#' @title A function to generate a dataframe of mixed cell-type origin simulated
#' samples using given reference atlas.
#' @param numberOfSamples The number of simulated samples to be generated in
#' the dataframe.
#' @param mixingVector Specify the cell origin proportions.If
#' numberOfSamples = 1, this can be a vector of length = number of cell types
#' in reference.Otherwise, this is a dataframe with rows for cell types (must
#' be equal to cell types in reference) and columns for samples. Cells contain
#' the proportion of the sample from the cell type. Use zeros for any unused
#' cell type. If this object is not given,will use random values for the
#' simulation.
#' @param reference A dataframe containing signatures of different cell types
#' used to generate the simulation. The first column should contain a unique
#' ID (e.g. CpG target ID) which can be used in deconvolution to match rows of
#' the reference to rows of the bulk. All subsequent columns are cell types.
#' Rows are units of the signature. Each cell contains the value for the cell
#' type and signature unit (e.g. methylation value at this CpG).
#' @importFrom methods is
#' @importFrom stats runif
#' @keywords simulation
#' @examples
#' data("HumanCellTypeMethAtlas")
#' bulk_mix50 <- simulateCellMix(50, reference = HumanCellTypeMethAtlas)
#'
#' bulk_mixVec <- simulateCellMix(1, mixingVector = c(
#'   0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#'   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#' ), reference = HumanCellTypeMethAtlas)
#' @return A list containing two data frames. simulated: A dataframe which
#' contains mixed cell-type origin simulated samples. The first column contains
#' a unique ID (used from reference) which can be used in deconvolution to match
#' rows of the reference to rows of the bulk.All subsequent columns are cell
#' types. Rows are units of signature (e.g. CpGs). Each cell contains the
#' value for the cell type and unit (e.g. methylation value at this CpG)
#' proportions: A dataframe with the cell proportions of the generated samples.
#' Each row is a sample. Columns are cell types.
#' @references Moss, J. et al.  (2018). Comprehensive human cell-type
#' methylation atlas reveals origins of circulating cell-free DNA in health
#' and disease. Nature communications, 9(1), 1-12.
#'   \url{https://doi.org/10.1038/s41467-018-07466-6}
#' @export

simulateCellMix <- function(numberOfSamples, mixingVector = NULL,
                            reference) {
  reference <- as.data.frame(reference)                            
  simulatedMixtureTable <- data.frame(matrix(
    ncol = numberOfSamples + 1,
    nrow = nrow(reference)
  ))
  ## setting up simulatedMixtureTable will hold simulation methylation values
  simulatedMixtureTable[, ] <- 0
  simulatedMixtureTable[, 1] <- reference[, 1]
  # copy IDs from reference atlas
  colnames(simulatedMixtureTable)[1] <- "IDs"

  proportionsTable <- data.frame(matrix(
    ncol = ncol(reference) - 1,
    nrow = numberOfSamples
  ))
  ## setting up proportionsTable will hold simulation origin cell proportions
  proportionsTable[, ] <- 0
  colnames(proportionsTable) <- colnames(reference[, -1])

  if (is.null(mixingVector)) {
    for (i in seq_len(numberOfSamples)) {
      n <- sample(seq(ncol(reference[-1])), 1)
      ## n is the randomly chosen number of cell origins
      picks <- sample(seq(ncol(reference[-1])), n, replace = FALSE)
      amts <- runif(n)
      ## amts is the randomly chosen proportions
      amts <- amts / sum(amts)
      ## normalize so the proportions add to 1

      for (c in seq_len(n)) {
        simulatedMixtureTable[, i + 1] <-
          (simulatedMixtureTable[, i + 1] +
            (amts[c] * reference[, picks[c] + 1]))
        # add the influence to the values
        proportionsTable[i, picks[c]] <- amts[c]
      }
      rownames(proportionsTable)[i] <- paste("Sample", i)
      colnames(simulatedMixtureTable)[i + 1] <- paste("Sample", i)
    }
  } else {
    if (is(mixingVector, "data.frame")) {
      if (ncol(mixingVector) != numberOfSamples) {
        stop("numberOfSamples must have the same number of rows with the
                mixingVector")
      }
      if (nrow(mixingVector) != (ncol(reference) - 1)) {
        stop("Number of cell types in mixingVector and reference must
                be equal (use zeros for unused cell types)")
      }
      if (typeof(unlist(mixingVector)) != "double") {
        stop("mixingVector should only contain numbers")
      } else {
        for (s in seq_len(ncol(mixingVector))) {
          # each row is a sample
          for (t in seq_len(nrow(mixingVector))) {
            # each column is a cell type
            if (mixingVector[t, s] > 0) {
              simulatedMixtureTable[, s + 1] <-
                (simulatedMixtureTable[, s + 1]
                + (mixingVector[t, s] * reference[, t + 1]))
              # add the influence to the values
              proportionsTable[s, t] <- mixingVector[t, s]
            }
          }
          rownames(proportionsTable)[s] <- paste("Sample", s)
          colnames(simulatedMixtureTable)[s + 1] <- paste("Sample", s)
        }
      }
    } else if (is(mixingVector, "numeric")) {
      if (1 != numberOfSamples) {
        stop("Only use a vector for mixingVector if numberOfSamples = 1,
            otherwise use dataframe with columns for samples and rows for cell
            type")
      }
      if (typeof(mixingVector) != "double") {
        stop("mixingVector should only contain numbers")
      }
      if (length(mixingVector) != (ncol(reference) - 1)) {
        stop("length of mixingVector must equal number of samples in
            reference (use zeros for unused cell types)")
      } else {
        for (t in seq_along(mixingVector)) {
          # mixing vector has number per cell type
          if (mixingVector[t] > 0) {
            simulatedMixtureTable[, 2] <-
              (simulatedMixtureTable[, 2] +
                (mixingVector[t] * reference[, t + 1]))
            # add the influence to the values
            proportionsTable[1, t] <- mixingVector[t]
          }
        }
        rownames(proportionsTable)[1] <- "Sample 1"
        colnames(simulatedMixtureTable)[2] <- "Sample 1"
      }
    } else {
      stop("mixingVector must be either data frame or numeric vector")
    }
  }
  return(list(
    simulated = simulatedMixtureTable,
    proportions = proportionsTable
  ))
}
