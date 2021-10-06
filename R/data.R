#' A dataset Illumina probe IDs of 400000 genomic loci (identified using the
#' “seqnames”, “ranges”, and “strand” values).
#' @usage data(probe_ids)
#' @format GRanges object with 400000 ranges and 1 metadata column:
#' \describe{
#'   \item{metadata}{contains seqnames ranges strand and ID}
#'   ...
#' }
#' @source \url{https://support.illumina.com/downloads/infinium-methylationepic-v1-0-product-files.html}
"probe_ids"

#' An example of WGBS data has been included as the GRanges object
#' This object contains the methylation values of 11715 genomic loci
#' (identified using the “seqnames”, “ranges”, and “strand” values) for a
#' single example sample. It can be used directly as WGBS_data in BSmeth2Probe.
#' @usage data(WGBS_GRanges)
#' @format GRanges object with 11715 ranges and 1 metadata column:
#' \describe{
#'   \item{metadata}{contains seqnames ranges strand and ID.
#'   seqinfo: 25 sequences from an unspecified genome; no seqlengths}
#'   ...
#' }
#' @source \url{https://www.encodeproject.org/data-standards/wgbs/}
"WGBS_GRanges"

#' The comprehensive human methylome reference atlas
#' @usage data("HumanCellTypeMethAtlas")
#' @format data.frame object with 6000 CpG loci and 25 human cell types column:
#' \describe{
#'   \item{cell type columns}{For each cell type and CpG locus, a methylation
#'   valuebetween 0 and 1 is provided. This value represents the fraction of
#'   methylated bases of the CpG locus.}
#'   \item{IDs}{CpG loci IDs for each cell type.}
#'   ...
#' }
#' @references Moss, J. et al.  (2018). Comprehensive human cell-type
#' methylation atlas reveals origins of circulating cell-free DNA in health
#' and disease. Nature communications, 9(1), 1-12.
#' @source \url{https://doi.org/10.1038/s41467-018-07466-6}
"HumanCellTypeMethAtlas"
