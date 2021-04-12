#' A function to construct a signature matrix
#' @param samples dataframe, has first column IDs, rest of columns are samples (must have column name as sample accession ID which should be found in sampleMeta), rows are units of signature (e.g. CpGs)
#' @param sampleMeta dataframe, must have Experiment_accession for accession ID of each sample, and Biosample_term_name for cell type of sample, rows are samples
#' @param atlas dataframe, the reference atlas to which new signatures can be added, if not present then a new reference atlas will be created using sample(s). Should be dataframe with column for each cell type, rows units of signature (e.g. CpGs)
#' @param variation_cutoff either a number between 0 to 1, or NULL. for multiple samples from the same cell type, ignore CpGs with variation > variation_cutoff with that cell type. defaults to NULL (i.e. no cutoff)
#' @importFrom magrittr %>%
#' @import data.table
#' @examples
#'  exampleSamples = simulateCellMix(1)[[1]]
#'  exampleReference = readRDS(system.file("reference_atlas_nodup.RDS", package = "deconvR"))
#'  exampleMeta = data.table( "Experiment_accession" = "example_sample",
#'                            "Biosample_term_name" = "example_cell_type")
#'  colnames(exampleSamples)[-1] =  c("example_sample")
#'  findSignatures(samples = exampleSamples, sampleMeta = exampleMeta)
#'  findSignatures(samples = exampleSamples, sampleMeta = exampleMeta, atlas = exampleReference)
#'
#' @return A dataframe extendedAtlas which contains all cell types in atlas (if given), and those in samples added by cell type, has first column "IDs", rest of columns are cell types, rows are have first cell with the ID (e.g. CpG ID) and then values of signature (e.g. methylation values)
#' @export

findSignatures = function(samples, sampleMeta, atlas=NULL, variation_cutoff = NULL) {

  if (!(is.null(variation_cutoff))) { # first just checking variation_cutoff is valid number if not null
    assertthat::assert_that(is.double(variation_cutoff), msg = "Variation cutoff should be a number")
    assertthat::assert_that(variation_cutoff >= 0, variation_cutoff <= 1, msg = "Variation cutoff should be between 0 and 1")
  }
  assertthat::assert_that(colnames(samples)[1] == "IDs", msg = "First column of samples should be IDs")
  for (sample_id in colnames(samples)[-1]) {
    assertthat::assert_that(sample_id %in% sampleMeta$Experiment_accession, msg = "All column names of samples (other than IDs column) should be present in sampleMeta Experiment_accession")
  }

  if(is.null(atlas)) {
    extendedAtlas = samples # if there's no atlas given then just work with samples
  }
  else {
    assertthat::assert_that(colnames(atlas)[1] == "IDs", msg = "First column of atlas should be IDs")
    for (celltype in colnames(atlas)[-1]) {                    # looking at all cell types in the atlas
        if (!(celltype %in% unlist(sampleMeta[,1]))) {
        sampleMeta  %>% tibble::add_row("Experiment_accession" = celltype, "Biosample_term_name" = celltype) # if the cell type isn't in metadata, just add a row saying it maps to itself
        }
      }
    extendedAtlas = data.table::merge.data.table(atlas,samples,by = "IDs") # merge atlas with samples
  }
  rowsToDelete = c()
  for (i in 1:NROW(sampleMeta)) {
    if (!(sampleMeta[i,1] %in% colnames(extendedAtlas))) {
      rowsToDelete = append(rowsToDelete, i)
    }
  }
  if (length(rowsToDelete) > 0) {
    sampleMeta = sampleMeta[-rowsToDelete,] # get rid of metadata entries not relevant to the atlas
  }

  #data.table::setDT(extendedAtlas)

  # units with variance across the entire atlas below 0.1% are excluded
  #rv = extendedAtlas[,matrixStats::rowVars(as.matrix(.SD)), .SDcols = names(extendedAtlas[1,!"IDs"])]
  #filteredExtendedAtlas = data.table::copy(extendedAtlas)
  #filteredExtendedAtlas = filteredExtendedAtlas[rv >= 0.001,]

  cat("CELL TYPES ADDED: \n")
  for (tissue in unique(sampleMeta$Biosample_term_name) ) {

      accession_ids = sampleMeta[Biosample_term_name %in% tissue, unique(Experiment_accession)]
      if (length(accession_ids) > 0) { # if no accession IDs map to this cell type then no need to continue
        tissueLabel = gsub(" ","_",tissue)
        if (!(is.null(variation_cutoff))) {
          # units with variance within tissue above variation_cutoff are excluded
          tissue_variance = extendedAtlas[,matrixStats::rowVars(as.matrix(.SD)), .SDcols = accession_ids]
          #tissue_variance = matrixStats::rowVars(as.matrix(extendedAtlas[,c(accession_ids)]))
          tissue_variance = nafill(tissue_variance, fill = 0) # variance NA (e.g. if only one sample), will be replaced with 0
          extendedAtlas = extendedAtlas[tissue_variance <= variation_cutoff,]
          }
        setDT(extendedAtlas)
        extendedAtlas = stats::na.omit(extendedAtlas) # units with missing values are excluded
        pooled_column = rowMeans(extendedAtlas[,.SD, .SDcols=accession_ids, drop = FALSE]) # pool samples for same tissue
        extendedAtlas[,(as.character(accession_ids)) := NULL ] # remove single samples
        extendedAtlas = cbind(extendedAtlas, pooled_column) # add new column pooling samples
        colnames(extendedAtlas)[NCOL(extendedAtlas)] = tissueLabel # pooled column named as Biosample_term_name
        cat(tissue,"\n")
      }
    }
  return(extendedAtlas)
}
