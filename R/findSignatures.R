#' A function to construct a signature matrix
#' @param samples dataframe, has first column "CpGs", rest of columns are samples (must have column name as sample accession ID which should be found in sampleMeta), rows are CpGs
#' @param sampleMeta dataframe, must have Experiment_accession for accession ID of each sample, and Biosample_term_name for cell type of sample, rows are samples
#' @param atlas dataframe, the reference atlas to which new signatures can be added, , if not present then a new reference atlas will be created using sample(s)
#' @keywords signaturehas first column "CpGs", rest of columns are cell types, rows are CpGs
#' @examples
#'  exampleAtlas_1 = readRDS(system.file("reference_atlas_nodup.RDS", package = "deconvR"))
#'  exampleMeta_1  = data.table("Experiment_accession" = as.character(1:25), "Biosample_term_name" = "example_cell_type")
#'  colnames(exampleAtlas_1)[2:26] =  as.character(1:25)
#'  findSignatures(samples = exampleAtlas_1, sampleMeta = exampleMeta_1)
#'
#'  exampleAtlas_2 = readRDS(system.file("reference_atlas_nodup.RDS", package = "deconvR"))
#'  exampleMeta_2  = data.table("Experiment_accession" = as.character(1:25), "Biosample_term_name" = colnames(exampleAtlas_2)[-1])
#'  exampleSamples_2 = exampleAtlas_2[,c(1,26)]
#'  exampleAtlas_2 = exampleAtlas_2[,c(1:25)]
#'  exampleMeta_2    = exampleMeta_2[25,]
#'  colnames(exampleSamples_2)[2]   = "25"
#'  findSignatures(samples = exampleSamples_2, sampleMeta = exampleMeta_2, atlas = exampleAtlas_2)
#'
#' @return A dataframe filteredExtendedMethAtlas which contains all cell types in atlas (if given), and those in samples added by cell type, has first column "CpGs", rest of columns are cell types, rows are CpGs, cells are methylation values
#' @export

findSignatures = function(samples, sampleMeta, atlas=NULL) {
  if(is.null(atlas)) {
    #extendedAtlas = data.table:merge.data.table(samples["CpGs"],samples,by = "CpGs")
    extendedAtlas = samples

  }
  else {
    extendedAtlas = data.table::merge.data.table(atlas,samples,by = "CpGs")
  }

  setDT(extendedAtlas)

  # CpGs with variance across the entire methylation atlas below 0.1% are excluded
  rv = extendedAtlas[, matrixStats::rowVars(as.matrix(.SD)), .SDcols = names(extendedAtlas[1,!"CpGs"])]
  filteredExtendedAtlas = data.table::copy(extendedAtlas)
  filteredExtendedAtlas = filteredExtendedAtlas[rv >= 0.001,]
  # CpGs with missing values are excluded
  filteredExtendedAtlas = stats::na.omit(filteredExtendedAtlas)
  cat("CELL TYPES ADDED: \n")
  for (tissue in unique(sampleMeta$Biosample_term_name) ) {
    cat(tissue,"\n")
    accession_ids = as.character(sampleMeta[Biosample_term_name %in% tissue, unique(Experiment_accession)])
    tissueLabel = gsub(" ","_",tissue)
    # http://brooksandrew.github.io/simpleblog/articles/advanced-data-table/#assign-a-column-with--named-with-a-character-object
    # pool samples for same tissue
    filteredExtendedAtlas[,(as.character(tissueLabel)) := rowMeans(.SD) , .SDcols = accession_ids]
    # remove single samples
    filteredExtendedAtlas[,(as.character(accession_ids)) := NULL ]
  }
  return(filteredExtendedAtlas)
}
