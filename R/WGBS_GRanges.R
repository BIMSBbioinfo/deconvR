#' An example of WGBS data has been included as the GRanges object
#' This object contains the methylation values of 11715 genomic loci
#' (identified using the “seqnames”, “ranges”, and “strand” values) for a single
#'  example sample. It can be used directly as WGBS_data in BSmeth2Probe.
#' @format GRanges object with 11715 ranges and 1 metadata column:
#' \describe{
#'   \item{metadata}{contains seqnames ranges strand and ID.
#'   seqinfo: 25 sequences from an unspecified genome; no seqlengths}
#'   ...
#' }
#' @source \url{https://www.encodeproject.org/data-standards/wgbs/}
"WGBS_GRanges"
