#' @title A function to map either multiple or single CpGs based on overlaps
#' @importFrom IRanges  mergeByOverlaps distance
#' @importFrom stats aggregate
#' @importFrom BiocGenerics as.data.frame
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
#' @param multipleMapping When searching for matches for probes not directly
#' covered in WGBS data, should WGBS CpGs which have already been mapped to
#' another probe still be considered? If TRUE, then yes. If FALSE, then no.
#' @param cutoff The maximum number of basepairs distance to consider for probes
#' which have not been directly covered in the WGBS data. Default value is 10.
#' @param overlaps_df A data frame containing the WGBS data and probe IDs
#' combined by overlaps
#' @return A dataframe with column for CpG IDs containing all of the maped CpGs 
#' based on overlaps
#' @noRd
mapByOverlaps <- function(WGBS_data, probe_id_locations, cutoff,
    multipleMapping) {
    overlaps_df <- mergeByOverlaps(WGBS_data, probe_id_locations)
    # mapping CpG locations to probe locations by overlap
    if (nrow(overlaps_df) == 0) {
        overlaps_df$distance <- numeric()
    }
    if (nrow(overlaps_df) > 0) {
        overlaps_df[, "distance"] <- NA
        # since there is no gap, say distance is NA
    }

    if (cutoff > 0) {
        # only need to do "nearlyOverlaps" if cutoff > 0
        allresults <- mapByNearOverlaps(
            WGBS_data, probe_id_locations, cutoff,
            multipleMapping, overlaps_df
        )
    }
    if (cutoff == 0) {
        ## if cutoff is 0, then all results are just exact overlaps
        allresults <- overlaps_df
    }
    ## cleaning up
    allresults <- allresults[, -c(1, NCOL(allresults) - 2, NCOL(allresults))]
    allresults <- allresults[c(ncol(allresults), seq_len(ncol(allresults) - 1))]
    ## if a probe was mapped to multiple CpGs,take the mean  value
    if (nrow(allresults) > 0) {
        allresults <- aggregate(
            x = allresults[, -1],
            by = list(ID = allresults[, 1]),
            FUN = mean
        )
        ## set column names to match WGBS data
        for (i in seq(2, ncol(allresults))) {
            colnames(allresults)[i] <-
                colnames(mcols(WGBS_data))[i - 1]
        }
    }
    ## the return value is a dataframe with column for CpG IDs
    colnames(allresults)[1] <- "IDs"
    allresults <- as.data.frame(allresults)
}
