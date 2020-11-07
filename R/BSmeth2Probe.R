#' A function to map WGBS methylation data to Illumina Probe IDs
#' @param probe_id_locations Either a dataframe or GRanges object containing probe IDs and their locations. If dataframe: must contain columns named "ID", "CHR", "Start", "End", and "Strand". If GRanges: should have locations ("seqnames", "ranges", "strand"), as well as metadata column "ID". Start and end locations should be 1-based coordinates. Note that any row with NA values will not be used. Example dataframe illumina_probes_hg38_GRanges.RDS in inst folder included in package.
#' @param WGBS_data Either a GRanges object or methylKit object (methylRaw, methylBase, methylRawDB, or methylBaseDB) of CpG locations and their methylation values. Contains locations ("seqnames", "ranges", "strand") and metadata column(s) of methylation values of sample(s) (i.e. one column per sample). These methylation values must be between 0 and 1.
#' @param cutoff The maximum number of basepairs distance to consider for probes which have not been directly covered in the WGBS data. Default value is 10.
#' @param multipleMapping When searching for matches for probes not directly covered in WGBS data, should WGBS CpGs which have already been mapped to another probe still be considered? If TRUE, then yes. If FALSE, then no. Default value is FALSE.
#' @keywords mapping
#' @examples
#' wgbs = readRDS(system.file("WGBS_GRanges.RDS", package = "deconvR"))
#' probe_ids = readRDS(system.file("illumina_probes_hg38_GRanges.RDS", package = "deconvR"))
#' BSmeth2Probe(probe_id_locations=probe_ids, WGBS_data=wgbs)
#' BSmeth2Probe(probe_id_locations=probe_ids, WGBS_data=wgbs[5:1000], cutoff = 0)
#' BSmeth2Probe(probe_id_locations=probe_ids, WGBS_data=wgbs, cutoff = 500, multipleMapping = TRUE)
#' BSmeth2Probe(probe_id_locations=probe_ids[100:200], WGBS_data=wgbs[1:10])
#' @return A dataframe with first column "CpGs" for CpG IDs, then 1 or more columns for methylation values of sample(s) (same number of samples as in WGBS_data)

#' ID for each probe which was mapped, and then methylation value(s) of the WGBS CpG to which it was matched (where either it overlapped or the gap was < cutoff).
#' If it matched to more than one CpG, the mean methylation value is taken.
#' @export

BSmeth2Probe = function(probe_id_locations, WGBS_data, cutoff = 10, multipleMapping = FALSE) {
  if (cutoff < 0) {stop("cutoff must be >= 0")}
  if (NROW(probe_id_locations) == 0) {stop("probe_id_locations should not be empty")}
  if (NROW(WGBS_data) == 0) {stop("WGBS_data should not be empty")}

  if ((class(probe_id_locations) != "data.frame") && (class(probe_id_locations) != "GRanges")) {
    stop("probe_id_locations must be either dataframe or GRanges object")
  }

  if (class(WGBS_data) != "GRanges" && class(WGBS_data) != "methylRawDB"&& class(WGBS_data) != "methylRaw" && class(WGBS_data) != "methylBaseDB"&& class(WGBS_data) != "methylBase") {
    stop("WGBS_data must be either GRanges object or methylKit object (methylRaw, methylBase, methylRawDB, or methylBaseDB")
  }
  if (class(probe_id_locations) == "data.frame") {
    if (any(is.null(probe_id_locations$CHR), is.null(probe_id_locations$Start), is.null(probe_id_locations$End), is.null(probe_id_locations$Strand), is.null(probe_id_locations$ID))) {
      stop("probe_id_locations must contain columns named ID, CHR, Start, End, and Strand.")
    }
    if (any(is.na(probe_id_locations[,c("CHR","Start","End","Strand","ID")]))) {
      print(paste("Dropping row containing NA: "+which(is.na(probe_id_locations))))
      probe_id_locations = tidyr::drop_na(probe_id_locations)
    }
    probe_id_locations = GenomicRanges::GRanges(seqnames=probe_id_locations[,"CHR"],  #if the probe_id_locations is a dataframe, turn it into GRanges before continuing
                                 ranges=IRanges::IRanges(start=probe_id_locations[,"Start"],
                                                end=probe_id_locations[,"End"]),
                                 strand=probe_id_locations[,"Strand"],
                                 ID = probe_id_locations[,"ID"])

  }
  if (class(probe_id_locations) =="GRanges") {
    if (is.null(probe_id_locations$ID)) {
      stop("probe_id_locations must have metadata column ID")
    }
  }
  if ((class(WGBS_data) == "methylBaseDB") || (class(WGBS_data) == "methylBase")) {
    percMeths = methylKit::percMethylation(WGBS_data)
    percMeths = percMeths / 100
    WGBS_data = methods::as(WGBS_data, "GRanges")
    GenomicRanges::mcols(WGBS_data) = as.data.frame(percMeths)
  }

  if (class(WGBS_data) != "GRanges") { #turn WGBS_data to GRanges object if it is not already one
    sampleName = methylKit::getSampleID(WGBS_data)
    WGBS_data = methods::as(WGBS_data, "GRanges")
    GenomicRanges::mcols(WGBS_data) = WGBS_data$numCs / WGBS_data$coverage
    colnames(GenomicRanges::mcols(WGBS_data)) = sampleName
  }
  if (class(WGBS_data) == "GRanges") {
    if (NCOL(S4Vectors::elementMetadata(WGBS_data)) == 0) {
      stop("WGBS_data should have at least one metadata column. No samples?")
    }
    if (any(BiocGenerics::lapply(as.data.frame(S4Vectors::elementMetadata(WGBS_data)), class)  != "numeric")) {
      stop("WGBS_data metadata columns should be methylation values of sample(s), between 0 and 1")
    }
    if ((any(as.data.frame(S4Vectors::elementMetadata(WGBS_data)) < 0)) || (any(as.data.frame(S4Vectors::elementMetadata(WGBS_data)) > 1))) {
      stop("WGBS_data metadata columns should be methylation values of sample(s), between 0 and 1")
    }
  }

  overlaps_df = IRanges::mergeByOverlaps(WGBS_data, probe_id_locations) #mapping CpG locations to probe locations by overlap
  #overlaps_df = cbind(overlaps_df, methylation=overlaps_df[,"numCs"]/overlaps_df[,"coverage"])  #methylation is fraction of coverage which read Cs
  if (nrow(overlaps_df) == 0) {overlaps_df$distance = numeric()}
  if (nrow(overlaps_df) > 0) {overlaps_df[,"distance"] <- NA} #since there is no gap, say distance is NA

  nearlyOverlaps_df = IRanges::mergeByOverlaps(WGBS_data, probe_id_locations, maxgap = cutoff)     #same mapping as first time, but now with cutoff gap allowed
  #nearlyOverlaps_df = cbind(nearlyOverlaps_df, methylation = nearlyOverlaps_df[,"numCs"]/nearlyOverlaps_df[,"coverage"]) #methylation is fraction of reads which were Cs
  nearlyOverlaps_df = subset(nearlyOverlaps_df, !(nearlyOverlaps_df$ID %in% overlaps_df$ID)) #discard probes already mapped in first round
  if (nrow(nearlyOverlaps_df) == 0) {nearlyOverlaps_df$distance = numeric()}
  if (nrow(nearlyOverlaps_df) > 0) {
    nearlyOverlaps_df = cbind(nearlyOverlaps_df, distance=IRanges::distance(nearlyOverlaps_df$WGBS_data, nearlyOverlaps_df$probe_id_locations)) #the distance of the gap between probe location and WGBS data location
    nearlyOverlaps_df = nearlyOverlaps_df[order(nearlyOverlaps_df$distance),] #order the df by distance so that when we delete duplicates, the duplicate with the largest gap is deleted

    if (!multipleMapping) {
      if (nrow(overlaps_df) > 0) {
        nearlyOverlaps_df = subset(nearlyOverlaps_df, !(is.element(data.table::transpose(BiocGenerics::as.data.frame(nearlyOverlaps_df$WGBS_data)), data.table::transpose(BiocGenerics::as.data.frame(overlaps_df$WGBS_data))))) #remove where CpG already mapped in first round if multipleMapping has been set to false
      }
      nearlyOverlaps_df = nearlyOverlaps_df[!duplicated(nearlyOverlaps_df$WGBS_data),]  #remove multiple mappings of CpG if multipleMapping has been set to false
    }
  }

  allresults = rbind(overlaps_df, nearlyOverlaps_df) #combine all of the mapping into one dataframe
  allresults = allresults[,-c(1,NCOL(allresults)-2, NCOL(allresults))]
  allresults = allresults[c(ncol(allresults), 1:ncol(allresults)-1)]
  if (nrow(allresults) > 0) {
    allresults = stats::aggregate(x = allresults[,-1], by = list(ID = allresults[,1]), FUN = mean) #if a probe was mapped to multiple CpGs, take the mean methylation value
    for (i in 2:ncol(allresults)) {
      colnames(allresults)[i] = colnames(GenomicRanges::mcols(WGBS_data)[i-1])
    }
  }
  colnames(allresults)[1] = "CpGs"
  if (nrow(allresults) == 0) {
    print("Result dataframe is empty. No matches could be found.")

        }

  as.data.frame(allresults) #the return value is a dataframe with column for CpG IDs, then 1 or more columns for methylation values of sample(s)

}