#' A function to map WGBS methylation data to Illumina Probe IDs
#' @param probe_id_locations Either a dataframe or GRanges object containing
#' probe IDs and their locations. If dataframe: must contain columns named "ID",
#' "seqnames", "Start", "End", and "Strand". If GRanges: should have locations
#' ("seqnames", "ranges", "strand"), as well as metadata column "ID". Start and
#' end locations should be 1-based coordinates. Note that any row with NA values
#' will not be used. Example dataframe illumina_probes_hg38_GRanges.RDS in inst
#' folder included in package.
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
#' @importFrom methods is isClass as
#' @importFrom tidyr drop_na
#' @importFrom stats aggregate
#' @importFrom methylKit percMethylation unite getSampleID
#' @importFrom GenomicRanges GRanges
#' @importFrom data.table nafill transpose
#' @importFrom IRanges IRanges mergeByOverlaps distance
#' @importFrom S4Vectors 'elementMetadata<-' runValue elementMetadata
#' @importFrom S4Vectors 'mcols<-' runValue mcols
#' @importFrom BiocGenerics  lapply as.data.frame
#' @keywords mapping
#' @examples
#' data("probe_ids")
#' data("WGBS_GRanges")
#' meth_probres <- BSmeth2Probe(
#'     probe_id_locations = probe_ids,
#'     WGBS_data = WGBS_GRanges
#' )
#' methp_cut <- BSmeth2Probe(
#'     probe_id_locations = probe_ids, WGBS_data = WGBS_GRanges[5:1000],
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
    if (NROW(probe_id_locations) == 0) {
        stop("probe_id_locations should not be empty")
    }
    if (NROW(WGBS_data) == 0) {
        stop("WGBS_data should not be empty")
    }

    if (!(class(probe_id_locations) %in% c(
        "data.frame", "GRanges"
    ))) {
        stop("Probe IDs must be either dataframe or GRanges objects.")
    }

    if (!(class(WGBS_data) %in% c(
        "GRanges", "methylRawList",
        "methylBaseDB", "methylBase"
    ))) {
        stop("WGBS_data must be either GRanges object or, a methylKit object,
            such as methylBase, methylRawList, methylBaseDB")
    }

    if (is.data.frame(probe_id_locations) ) {
        names(probe_id_locations) <- vapply(names(probe_id_locations), tolower,
            FUN.VALUE = character(1)
        )
        if (any(is.null(c("seqnames", "end", "start", "strand", "id") %in%
            names(probe_id_locations))) ) {
            stop("probe_id_locations must contain columns named ID, Seqnames,
                 Start, End, and Strand.")
        }
        if (anyNA(probe_id_locations) ) {
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
    if (is(probe_id_locations, "GRanges") ) {
        if (is.null(probe_id_locations$ID) ) {
            stop("probe_id_locations must have metadata column ID")
        }
    }
    if ((is(WGBS_data, "methylBaseDB") ) ||
        (is(WGBS_data, "methylBase") )) {
        pm <- percMethylation(WGBS_data)
        pm <- pm / 100
        WGBS_data <- as(WGBS_data, "GRanges")
        mcols(WGBS_data) <- as.data.frame(pm)
    }
    if (is(WGBS_data, "methylRawList") ) {
        ## Convert methylRawList to methylbase object
        WGBS_data <- unite(WGBS_data, destrand = FALSE)
        pm <- percMethylation(WGBS_data)
        pm <- pm / 100
        WGBS_data <- as(WGBS_data, "GRanges")
        mcols(WGBS_data) <- as.data.frame(pm)
    }
    ## turn WGBS_data to GRanges object if it is not already one
    if (isClass(WGBS_data, Class = "GRanges") != TRUE) {
        sampleName <- getSampleID(WGBS_data)
        WGBS_data <- as(WGBS_data, "GRanges")
        mcols(WGBS_data) <- WGBS_data$numCs / WGBS_data$coverage
        colnames(mcols(WGBS_data)) <- sampleName
    }
    if (isClass(WGBS_data, Class = "GRanges") ) {
        column_names <- colnames(mcols(WGBS_data))
        elementMetadata(WGBS_data) <-
            nafill(as.data.frame(elementMetadata(WGBS_data)),
                fill = 0
            )
        colnames(mcols(WGBS_data)) <- column_names
        if (NCOL(elementMetadata(WGBS_data)) == 0) {
            stop("WGBS_data must have at least one metadata column.")
        }
        if (any(lapply(
            as.data.frame(elementMetadata(WGBS_data)),
            class
        ) != "numeric")) {
            stop("The metadata columns of WGBS_data must contain methylation
                 values of sample(s) between 0 and 1")
        }
        if ((any(as.data.frame(elementMetadata(WGBS_data)) < 0)) ||
            (any(as.data.frame(elementMetadata(WGBS_data)) > 1))) {
            stop("The metadata columns of WGBS_data must contain methylation
                values of sample(s) between 0 and 1")
        }
    }

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
            ## the distance of the gap between probe  and WGBS data location
            nearolaps_df <- nearolaps_df[order(nearolaps_df$distance), ]
            # order the df by distance so that when we delete duplicates
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
                ## remove multiple mappings of CpG if multipleMapping is false
                nearolaps_df <-
                    nearolaps_df[!duplicated(nearolaps_df$WGBS_data), ]
            }
        }
        allresults <- rbind(overlaps_df, nearolaps_df)
        ## combine all of the mapping into one dataframe
    }
    if (cutoff == 0) {
        ## if cutoff is 0, then all results are just exact overlaps
        allresults <- overlaps_df
    }

    allresults <- allresults[, -c(1, NCOL(allresults) - 2, NCOL(allresults))]
    ## these are just cleaning up allresults a bit to make aggregation easier
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
    colnames(allresults)[1] <- "IDs"
    if (nrow(allresults) == 0) {
        message("Result dataframe is empty. No matches could be found.")
    }
    ## the return value is a dataframe with column for CpG IDs
    return(as.data.frame(allresults))
}
