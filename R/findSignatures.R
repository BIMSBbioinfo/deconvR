#' A function to construct a signature matrix
#' @param samples dataframe, has first column for IDs, rest of columns are samples (must have column name as sample accession ID which should be found in sampleMeta), rows are units of signature (e.g. CpGs)
#' @param sampleMeta dataframe, must have Experiment_accession for accession ID of each sample, and Biosample_term_name for cell type of sample, rows are samples
#' @param atlas dataframe, the reference atlas to which new signatures can be added, if not present then a new reference atlas will be created using sample(s). Should be dataframe with column for each cell type, rows units of signature (e.g. CpGs)
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
#' @return A dataframe filteredExtendedAtlas which contains all cell types in atlas (if given), and those in samples added by cell type, has first column "IDs", rest of columns are cell types, rows are have first cell with the ID (e.g. CpG ID) and then values of signature (e.g. methylation values)
#' @export

findSignatures = function(samples, sampleMeta, atlas=NULL) {
  if(is.null(atlas)) {
    extendedAtlas = samples # if there's no atlas given then just work with samples
  }
  else {
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

  data.table::setDT(extendedAtlas)

  # units with variance across the entire atlas below 0.1% are excluded
  rv = extendedAtlas[, matrixStats::rowVars(as.matrix(.SD)), .SDcols = names(extendedAtlas[1,!"IDs"])]
  filteredExtendedAtlas = data.table::copy(extendedAtlas)
  filteredExtendedAtlas = filteredExtendedAtlas[rv >= 0.001,]
  # units with missing values are excluded
  filteredExtendedAtlas = stats::na.omit(filteredExtendedAtlas)
  cat("CELL TYPES ADDED: \n")
  for (tissue in unique(sampleMeta$Biosample_term_name) ) {
      accession_ids = as.character(sampleMeta[sampleMeta$Biosample_term_name %in% tissue, unique(Experiment_accession)])
      if (length(accession_ids) > 0) { # if no accession IDs map to this cell type then no need to continue
        tissueLabel = gsub(" ","_",tissue)
        # http://brooksandrew.github.io/simpleblog/articles/advanced-data-table/#assign-a-column-with--named-with-a-character-object
        # pool samples for same tissue
        filteredExtendedAtlas[,(as.character(tissueLabel)) := rowMeans(.SD) , .SDcols = accession_ids]
        # remove single samples
        filteredExtendedAtlas[,(as.character(accession_ids)) := NULL ]
        cat(tissue,"\n")
      }
    }
  return(filteredExtendedAtlas)
}
