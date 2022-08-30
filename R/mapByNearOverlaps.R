#' @title A function to map either multiple or single CpGs based on near
#' overlaps if cutoff > 0
#' @importFrom IRanges  mergeByOverlaps distance
#' @importFrom data.table nafill transpose
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
#' @keywords internal
#' @return A dataframe containing all of the maped CpGs based on near overlaps
#' @noRd
mapByNearOverlaps <- function(WGBS_data, probe_id_locations, cutoff,
                              multipleMapping, overlaps_df) {
  # only need to do "nearlyOverlaps" if cutoff > 0
  nearolaps_df <- mergeByOverlaps(WGBS_data,
    probe_id_locations,
    maxgap = cutoff
  )
  ## same mapping as first time, but now with cutoff gap allowed
  nearolaps_df <- subset(
    nearolaps_df,
    !(nearolaps_df$ID %in% overlaps_df$ID)
  )
  ## discard probes already mapped in first round
  if (nrow(nearolaps_df) == 0) {
    nearolaps_df$distance <- numeric()
  }
  if (nrow(nearolaps_df) > 0) {
    nearolaps_df <- cbind(nearolaps_df,
      distance = distance(
        nearolaps_df$WGBS_data,
        nearolaps_df$probe_id_locations
      )
    )
    ## the distance of the gap between probe and WGBS data location
    nearolaps_df <- nearolaps_df[order(nearolaps_df$distance), ]
    # the duplicate with the largest gap is deleted
    if (!multipleMapping) {
      if (nrow(overlaps_df) > 0) {
        ## remove where CpG already mapped in first round
        ## if multipleMapping  has been set to false
        nearolaps_df <- subset(nearolaps_df, !(is.element(
          transpose(as.data.frame(
            nearolaps_df$WGBS_data
          )),
          transpose(as.data.frame(
            overlaps_df$WGBS_data
          ))
        )))
      }
      nearolaps_df <-
        nearolaps_df[!duplicated(nearolaps_df$WGBS_data), ]
    }
  }
  allresults <- rbind(overlaps_df, nearolaps_df)
  ## combine all of the mapping into one dataframe
}
