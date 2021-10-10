#' @title A function to map WGBS methylation data to Illumina Probe IDs
#' @param probe_id_locations Either a dataframe or GRanges object containing
#' probe IDs and their locations. If dataframe: must contain columns named "ID",
#' "seqnames", "Start", "End", and "Strand". If GRanges: should have locations
#' ("seqnames", "ranges", "strand"), as well as metadata column "ID". Start and
#' end locations should be 1-based coordinates. Note that any row with NA values
#' will not be used.
#' @param WGBS_data Either a GRanges object or methylKit object (methylRaw,
#' methylBase, methylRawDB, or methylBaseDB) of CpG locations and their
#' methylation values. Contains locations ("seqnames", "ranges", "strand") and
#' metadata column(s) of methylation values of sample(s) (i.e. one column per
#' sample). These methylation values must be between 0 and 1.
#' @param cutoff The maximum number of basepairs distance to consider for probes
#' which have not been directly covered in the WGBS data. Default value is 10.
#' @param multipleMapping When searching for matches for probes not directly
#' covered in WGBS data, should WGBS CpGs which have already been mapped to
#' another probe still be considered? If TRUE, then yes. If FALSE, then no.
#' Default value is FALSE.
#' @importFrom methods is as
#' @importFrom tidyr drop_na
#' @importFrom GenomicRanges GRanges
#' @importFrom data.table nafill
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors 'elementMetadata<-' runValue elementMetadata
#' @importFrom S4Vectors 'mcols<-' runValue mcols
#' @keywords mapping
#' @examples
#' data("IlluminaMethEpicB5ProbeIDs")
#' WGBS_GRanges <- readRDS(system.file("extdata", "WGBS_GRanges.RDS",
#'     package = "deconvR"
#' ))
#' meth_probres <- BSmeth2Probe(
#'     probe_id_locations = IlluminaMethEpicB5ProbeIDs,
#'     WGBS_data = WGBS_GRanges
#' )
#' methp_cut <- BSmeth2Probe(
#'     probe_id_locations = IlluminaMethEpicB5ProbeIDs,
#'     WGBS_data = WGBS_GRanges[5:1000],
#'     cutoff = 2
#' )
#' @return A dataframe with first column "IDs" for CpG IDs, then 1 or more
#' columns for methylation values of sample(s) (same number of samples as in
#' WGBS_data)
#' ID for each probe which was mapped, and then methylation value(s) of the
#' WGBS CpG to which it was matched (where either it overlapped or the gap was
#' < cutoff). If it matched to more than one CpG, the mean methylation value
#' is taken.
#' @export

BSmeth2Probe <- function(probe_id_locations, WGBS_data, cutoff = 10,
    multipleMapping = FALSE) {
    if (cutoff < 0) {
        stop("cutoff must be >= 0")
    }
    if (NROW(probe_id_locations) == 0 || NROW(WGBS_data) == 0) {
        stop("Probe ID locations or WGBS data must not be empty")
    }

    if (!(class(probe_id_locations) %in% c(
        "data.frame", "GRanges"
    ))) {
        stop("Probe IDs must be either dataframe or GRanges objects.")
    }

    if (!(class(WGBS_data) %in% c(
        "GRanges", "methylRawList", "methylRaw", "methylRawListDB",
        "methylBaseDB", "methylBase"
    ))) {
        stop("WGBS_data must be either GRanges object or, a methylKit object,
        such as methylBase, methylRawList, methylBaseDB")
    }

    if (is.data.frame(probe_id_locations)) {
        names(probe_id_locations) <- vapply(names(probe_id_locations), tolower,
            FUN.VALUE = character(1)
        )
        if (any(is.null(c("seqnames", "end", "start", "strand", "id") %in%
            names(probe_id_locations)))) {
            stop("probe_id_locations must contain columns named ID, Seqnames,
            Start, End, and Strand.")
        }
        if (anyNA(probe_id_locations)) {
            message("Dropping row containing NA: " +
                which(is.na(probe_id_locations)))
            probe_id_locations <- drop_na(probe_id_locations)
        }
        probe_id_locations <- GRanges(
            seqnames = probe_id_locations[, "seqnames"],
            ## turn it into GRanges before continuing
            ranges = IRanges(
                start = probe_id_locations[, "start"],
                end = probe_id_locations[, "end"]
            ),
            strand = probe_id_locations[, "strand"],
            ID = probe_id_locations[, "id"]
        )
    }
    if (is(probe_id_locations, "GRanges")) {
        if (is.null(probe_id_locations$ID)) {
            stop("probe_id_locations must have metadata column ID")
        }
    }
    # Convert methylKit objects to GRanges
    if (!is(WGBS_data, "GRanges")) {
        WGBS_data <- convertMe(WGBS_data)
    }
    if (is(WGBS_data, "GRanges")) {
        column_names <- colnames(mcols(WGBS_data))
        elementMetadata(WGBS_data) <-
            nafill(as.data.frame(elementMetadata(WGBS_data)),
                fill = 0
            )
        colnames(mcols(WGBS_data)) <- column_names
        if (NCOL(elementMetadata(WGBS_data)) == 0) {
            stop("WGBS_data must have at least one metadata column.")
        }
        if ((any(as.data.frame(elementMetadata(WGBS_data)) < 0)) ||
            (any(as.data.frame(elementMetadata(WGBS_data)) > 1))) {
            stop("The metadata columns of WGBS_data must contain methylation
            values of sample(s) between 0 and 1")
        }
    }

    allresults <- mapByOverlaps(
        WGBS_data, probe_id_locations, cutoff,
        multipleMapping
    )
    if (nrow(allresults) == 0) {
        message("BSmeth2Probe couldn't find any match between the probe IDs and
        WGBS data you have provided. Please make sure that you're using an
        appropriate probe ID locations with your WGBS data.")
    }
    return(allresults)
}
