#' A function to generate a dataframe of mixed cell-type origin simulated samples using given reference atlas.
#' @param numberOfSamples The number of simulated samples to be generated in the dataframe.
#' @param reference A reference dataframe with which to generate the simulations. Should be dataframe with first column "CpGs" for Illumina IDs, and rest of columns = cell types. Should not have duplicate Illumina IDs. If not given, will default to ref_atlas included in package
#' @keywords simulation
#' @examples
#' makeBigTable(50)
#' makeBigTable(numberOfSamples=100,
#'   reference=read.csv(system.file("reference_atlas_nodup.csv",package = "deconvR")))
#' @return A dataframe "bigExampleTable" which contains mixed cell-type origin simulated samples. First column "CpGs" for Illumina IDs, and rest of columns = cell types.
#' @export

makeBigTable = function(numberOfSamples, reference=utils::read.csv(system.file("reference_atlas_nodup.csv", package = "deconvR"))) {
  bigExampleTable = data.frame(matrix(ncol=numberOfSamples+1,nrow=nrow(reference)))

  colnames(bigExampleTable) = "_"
  colnames(bigExampleTable)[1] = "CpGs"
  bigExampleTable[,] = 0
  bigExampleTable[,1] = reference[,1] #copy CpGs from reference atlas

for (i in 1:numberOfSamples) {
  n = sample(1:25, 1)  # n is the randomly chosen number of cell origins that will make up this sample
  amts = stats::runif(n)      # amts is the randomly chosen proportions of the different cell origins
  amts = amts/sum(amts) # normalize so the proportions add to 1

  for (samp in amts) {
    pick = sample(1:25, 1) # randomly pick the cell type origin
    if (!grepl(colnames(reference)[pick+1], colnames(bigExampleTable)[i+1], fixed = TRUE)) { # make sure we're not adding the same cell type twice
      bigExampleTable[,i+1] = (bigExampleTable[,i+1] + (samp * reference[,pick+1])) # add the influence to the values
      if (is.na(colnames(bigExampleTable)[i+1])){
        colnames(bigExampleTable)[i+1] = paste(colnames(reference)[pick+1],samp,"_", sep="") # naming the column to say the cellular type and proportion
      }
      else {
        colnames(bigExampleTable)[i+1] = paste(colnames(bigExampleTable)[i+1], colnames(reference)[pick+1],samp,"_", sep="")
      }
    }
  }
  colnames(bigExampleTable)[i+1] = substr(colnames(bigExampleTable)[i+1],1,nchar(colnames(bigExampleTable)[i+1])-1)

}
  return(bigExampleTable)
}
