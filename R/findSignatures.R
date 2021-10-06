#' @title A function to construct a signature matrix
#' @param samples dataframe, has first column IDs, rest of columns are samples
#' (must have column name as sample accession ID which should be found in
#' sampleMeta), rows are units of signature (e.g. CpGs)
#' @param sampleMeta dataframe, must have first column for accession ID of each
#' sample, and second column for cell type of sample, rows are samples
#' @param atlas dataframe, the reference atlas to which new signatures can be
#' added, if not present then a new reference atlas will be created using
#' sample(s). Should be dataframe with column for each cell type, rows units of
#' signature (e.g. CpGs)
#' @param variation_cutoff either a number between 0 to 1, or NULL.For multiple
#' samples from the same cell type, ignore CpGs with variation >
#' variation_cutoff with that cell type. defaults to NULL (i.e. no cutoff)
#' @importFrom magrittr %>%
#' @importFrom data.table  merge.data.table .SD setDT :=
#' @importFrom assertthat assert_that
#' @importFrom matrixStats rowVars
#' @importFrom stats na.omit
#' @examples
#' data("HumanCellTypeMethAtlas")
#' exampleSamples <- simulateCellMix(1, reference = HumanCellTypeMethAtlas)[[1]]
#' exampleMeta <- data.table(
#'     "Experiment_accession" = "example_sample",
#'     "Biosample_term_name" = "example_cell_type"
#' )
#' colnames(exampleSamples)[-1] <- c("example_sample")
#' signatures <- findSignatures(
#'     samples = exampleSamples,
#'     sampleMeta = exampleMeta
#' )
#' signatures <- findSignatures(
#'     samples = exampleSamples, sampleMeta = exampleMeta,
#'     atlas = HumanCellTypeMethAtlas
#' )
#' @return A dataframe extendedAtlas which contains all cell types in atlas
#' (if given), and those in samples added by cell type, has first column "IDs",
#' rest of columns are cell types, rows are have first cell with the ID
#' (e.g. CpG ID) and then values of signature (e.g. methylation values)
#' @export

findSignatures <- function(samples, sampleMeta, atlas = NULL,
    variation_cutoff = NULL) {
    if (!(is.null(variation_cutoff))) {
        ## first just checking variation_cutoff is valid number if not null
        assert_that(is.double(variation_cutoff),
            msg = "Variation cutoff should be a number"
        )
        assert_that(variation_cutoff >= 0,
            variation_cutoff <= 1,
            msg = "Variation cutoff should be between 0 and 1"
        )
    }
    assert_that(colnames(samples)[1] == "IDs",
        msg = "First column of samples should be IDs"
    )
    for (sample_id in colnames(samples)[-1]) {
        assert_that(sample_id %in% unlist(sampleMeta[, 1]),
            msg = "All column names of samples
                            (other than IDs column) should be present in
                            sampleMeta Experiment_accession"
        )
    }

    if (is.null(atlas)) {
        extendedAtlas <- samples
        ## if there's no atlas given then just work with samples
    } else {
        assert_that(colnames(atlas)[1] == "IDs",
            msg = "First column of atlas should be IDs"
        )
        for (celltype in colnames(atlas)[-1]) {
            if (!(celltype %in% unlist(sampleMeta[, 1]))) {
                sampleMeta <- rbind(sampleMeta, list(celltype, celltype))
                ## if the cell type isn't in metadata,
                # just add a row saying it maps to itself
            }
        }
        extendedAtlas <- merge.data.table(atlas, samples,
            by = "IDs", sort = FALSE
        )
        # merge atlas with samples
    }
    rowsToDelete <- c()
    for (i in seq_len(NROW(sampleMeta))) {
        if (!(sampleMeta[i, 1] %in% colnames(extendedAtlas))) {
            rowsToDelete <- append(rowsToDelete, i)
        }
    }
    if (length(rowsToDelete) > 0) {
        sampleMeta <- sampleMeta[-rowsToDelete, ]
        # get rid of metadata entries not relevant to the atlas
    }

    cat("CELL TYPES IN EXTENDED ATLAS: \n")
    for (tissue in unlist(unique(sampleMeta[, 2]))) {
        accession_ids <- sampleMeta[sampleMeta[[2]] %in% tissue, 1]
        if (length(accession_ids) > 0) {
            # if no accession IDs map to this cell type
            tissueLabel <- gsub(" ", "_", tissue)
            # units with variance within tissue above variation_cutoff excluded
            if (!(is.null(variation_cutoff))) {
                tissue_variance <-
                    extendedAtlas[, rowVars(as.matrix(.SD)),
                        .SDcols = unlist(accession_ids)
                    ]

                tissue_variance <- nafill(tissue_variance, fill = 0)
                # variance NA (e.g. if only one sample), will be replaced with 0
                extendedAtlas <- extendedAtlas[tissue_variance <=
                    variation_cutoff, ]
            }
            setDT(extendedAtlas)
            extendedAtlas <- na.omit(extendedAtlas)
            # units with missing values are excluded
            pooled_column <- rowMeans(extendedAtlas[, .SD,
                .SDcols = unlist(accession_ids),
                drop = FALSE
            ])
            # pool samples for same tissue
            extendedAtlas[, as.character(unlist(accession_ids)) := NULL]
            # remove single samples
            extendedAtlas <- cbind(extendedAtlas, pooled_column)
            # add new column pooling samples
            colnames(extendedAtlas)[NCOL(extendedAtlas)] <- tissueLabel
            # pooled column named as Biosample_term_name
            cat(tissue, "\n")
        }
    }
    return(extendedAtlas)
}
