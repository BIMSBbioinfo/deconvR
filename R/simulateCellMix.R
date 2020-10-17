#' A function to generate a dataframe of mixed cell-type origin simulated samples using given reference atlas.
#' @param numberOfSamples The number of simulated samples to be generated in the dataframe.
#' @param mixingVector Specify the cell origin proportions. If numberOfSamples = 1, this can be a vector of length = number of cell types in reference.
#' Otherwise, this is a dataframe with rows for cell types (must be equal to cell types in reference) and columns for samples.
#' Cells contain the proportion of the sample from the cell type. Use zeros for any unused cell type. If this object is not given,
#' will use random values for the simulation.
#' @param reference A dataframe containing CpG signatures of different cell types used to generate the simulation. The first column should contain a unique ID
#' (e.g. target ID) which can be used in deconvolution to match rows of the reference to rows of the bulk. All subsequent columns are cell types.
#' Rows are CpGs. Each cell contains the methylation value for the cell type at the CpG location. If not given, defaults to a reference atlas which is included in this
#' package (see deconvR/inst/reference_atlas_nodup.RDS). This atlas contains 25 cell types, 6105 CpGs (and their target IDs), and therefore 152,625 methylation values.
#' @keywords simulation
#' @examples
#' simulateCellMix(50)
#' simulateCellMix(numberOfSamples=100,
#'   reference=readRDS(system.file("reference_atlas_nodup.RDS", package = "deconvR")))
#' simulateCellMix(1,mixingVector = c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
#' simulateCellMix(1,data.frame( c(0,0,0,0,0,0,0,0.5,0,0,0,0,0,0,0,0,0,0,0,0.5,0,0,0,0,0)))
#' simulateCellMix(2,data.frame( c(0,0,0,0,0,0,0,0.5,0,0,0,0,0,0,0,0,0,0,0,0.5,0,0,0,0,0),
#'   c(0.1,0,0,0,0,0,0.5,0,0,0.2,0,0,0,0,0,0,0.2,0,0,0,0,0,0,0,0)))
#' @return A list containing two data frames.
#' First: A dataframe which contains mixed cell-type origin simulated samples. The first column contains a unique ID (used from reference) which can be used in
#' deconvolution to match rows of the reference to rows of the bulk. All subsequent columns are cell types. Rows are CpGs. Each cell contains the methylation
#' value for the cell type at the CpG location.
#' Second: A dataframe with the cell proportions of the generated samples. Each row is a sample. Columns are cell types.
#' @export

simulateCellMix = function(numberOfSamples, mixingVector=NULL, reference=readRDS(system.file("reference_atlas_nodup.RDS", package = "deconvR"))) {
  simulatedMixtureTable = data.frame(matrix(ncol=numberOfSamples+1,nrow=nrow(reference)))
#  colnames(simulatedMixtureTable) = "_"
  colnames(simulatedMixtureTable)[1] = "CpGs"
  simulatedMixtureTable[,] = 0
  simulatedMixtureTable[,1] = reference[,1] #copy CpGs from reference atlas

  proportionsTable = data.frame(matrix(ncol=ncol(reference)-1, nrow=numberOfSamples))
  proportionsTable[,] = 0
  colnames(proportionsTable) = colnames(reference[,-1])

  if (is.null(mixingVector)) {
    for (i in 1:numberOfSamples)  {
      n = sample(1:25, 1)  # n is the randomly chosen number of cell origins that will make up this sample
      picks = sample(1:25, n, replace = FALSE)
      amts = stats::runif(n)      # amts is the randomly chosen proportions of the different cell origins
      sumAmts=sum(amts)
      amts = amts/sumAmts # normalize so the proportions add to 1

      for (c in 1:n) {
          simulatedMixtureTable[,i+1] = (simulatedMixtureTable[,i+1] + (amts[c] * reference[,picks[c]+1])) # add the influence to the values
          proportionsTable[i, picks[c]] = amts[c]
      }
      rownames(proportionsTable)[i] = paste("Sample ", i)
      colnames(simulatedMixtureTable)[i+1] = paste("Sample ", i)

    }
  }
  else {
    if (class(mixingVector) == "data.frame") {
      if (ncol(mixingVector) != numberOfSamples) {stop("numberOfSamples should equal number of rows in mixingVector")}
      if (nrow(mixingVector) != (ncol(reference) - 1 )) {stop("number of cell types in mixingVector and reference should be equal (use zeros for unused cell types)")}
      if (typeof(unlist(mixingVector)) != "double") {stop("mixingVector should only contain numbers")}
      else {
        for (s in 1:ncol(mixingVector)) {  #each row is a sample
          for (t in 1:nrow(mixingVector))  { #each column is a cell type
            if (mixingVector[t,s] > 0) {
              simulatedMixtureTable[,s+1] = (simulatedMixtureTable[,s+1] + (mixingVector[t,s] * reference[,t+1])) # add the influence to the values
              proportionsTable[s, t] = mixingVector[t,s]
            }
          }
          rownames(proportionsTable)[s] = paste("Sample ", s)
          colnames(simulatedMixtureTable)[s+1] = paste("Sample ", s)

        }
      }
    }
    else if (class(mixingVector) == "numeric") {
      if (1 != numberOfSamples) {stop("you may only use a vector for mixingVector if numberOfSamples = 1, otherwise use dataframe with columns for samples and rows for cell type")}
      if (typeof(mixingVector) != "double") {stop("mixingVector should only contain numbers")}
      if (length(mixingVector) != (ncol(reference) -1)) {stop("length of mixingVector must equal number of samples in reference (use zeros for unused cell types)")}
      else {
        for (t in 1:length(mixingVector))  { #mixing vector has number per cell type
          if (mixingVector[t] > 0) {
            simulatedMixtureTable[,2] = (simulatedMixtureTable[,2] + (mixingVector[t] * reference[,t+1])) # add the influence to the values
            proportionsTable[1, t] = mixingVector[t]
          }
        }
        rownames(proportionsTable)[1] = "Sample 1"
        colnames(simulatedMixtureTable)[1] = "Sample 1"

      }
    }
    else {stop("mixingVector must be either data frame or numeric vector")}
  }
  return(list(simulatedMixtureTable, proportionsTable))
}









