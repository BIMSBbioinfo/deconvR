#' A function to generate a dataframe of mixed cell-type origin simulated samples using given reference atlas.
#' @param numberOfSamples The number of simulated samples to be generated in the dataframe.
#' @param reference A reference dataframe with which to generate the simulations. Should be dataframe with first column "CpGs" for Illumina IDs, and rest of columns = cell types. Should not have duplicate Illumina IDs.
#' @keywords simulation
#' @return A dataframe "bigExampleTable" which contains mixed cell-type origin simulated samples. First column "CpGs" for Illumina IDs, and rest of columns = cell types.
#' @export

#numberOfSamples = 500

makeBigTable = function(numberOfSamples, reference) {
  bigExampleTable = data.frame(matrix(ncol=numberOfSamples+1,nrow=nrow(reference)))

  colnames(bigExampleTable) = "_"
  colnames(bigExampleTable)[1] = "CpGs"
  bigExampleTable[,] = 0
  bigExampleTable[,1] = reference[,1]

for (i in 1:numberOfSamples) {
  n = sample(1:25, 1)
  amts = runif(n)
  amts = amts/sum(amts)

  for (samp in amts) {
    pick = sample(1:25, 1)
    bigExampleTable[,i+1] = (bigExampleTable[,i+1] + (samp * reference[,pick+1]))

    if (is.na(colnames(bigExampleTable)[i+1])){
      colnames(bigExampleTable)[i+1] = paste(colnames(reference)[pick+1],samp,"_", sep="")
    }
    else {
      colnames(bigExampleTable)[i+1] = paste(colnames(bigExampleTable)[i+1], colnames(reference)[pick+1],samp,"_", sep="")
    }
  }
  colnames(bigExampleTable)[i+1] = substr(colnames(bigExampleTable)[i+1],1,nchar(colnames(bigExampleTable)[i+1])-1)

}
  return(bigExampleTable)
}

# makeBigTable(500, reference_atlas_nodup)
