#' @title A function to calculate the percent of methylation and convert
#' methylkit objects to GRanges objects. Note that it does not support
#' methylDiff objects.
#' @param WGBS_data WGBS data in one of the following methylKit objects:
#' methylRawList, methylRaw, methylRawListDB,methylBaseDB, methylBase.
#' @importFrom methods is as
#' @importFrom dplyr select starts_with
#' @importFrom methylKit unite getSampleID
#' @importFrom GenomicRanges  makeGRangesFromDataFrame
#' @examples
#' library(methylKit)
#' data("methylKit") # Load the data
#' basegr <- convertMe(methylBase.obj)
#' rawlistgr <- convertMe(methylRawList.obj)
#' @return a GRanges object containing the methylation percent of samples
#' @keywords internal
#' @noRd
convertMe <- function(WGBS_data) {
  # Set sample IDs to use later
  sampleName <- getSampleID(WGBS_data)
  if (is(WGBS_data, "methylRawList") ||
    is(WGBS_data, "methylRawListDB")) {
    WGBS_data <- unite(WGBS_data, destrand = FALSE)
  }
  WGBS_data <- as(WGBS_data, "GRanges")
  WGBS_data <- as.data.frame(WGBS_data)
  # Calculate the percent of methylation
  percent_meth <- select(WGBS_data, starts_with("numCs")) /
    select(WGBS_data, starts_with("coverage"))
  WGBS_data <- makeGRangesFromDataFrame(WGBS_data,
    keep.extra.columns = TRUE
  )
  mcols(WGBS_data) <- percent_meth
  # set sample name
  colnames(mcols(WGBS_data)) <- sampleName
  return(WGBS_data)
}
