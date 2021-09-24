context("BSmeth2Probe")
library(deconvR)

test_that("BSmeth2Probe", {
    wgbs <- readRDS(system.file("WGBS_GRanges.RDS", package = "deconvR"))
    probe_ids <- readRDS(system.file("illumina_probes_hg38_GRanges.RDS", package = "deconvR"))
    probe_ids_df <- GenomicRanges::as.data.frame(probe_ids)

    results <- BSmeth2Probe(probe_id_locations = probe_ids, WGBS_data = wgbs)

    expect_equal(class(results), "data.frame")
    expect_equal(NCOL(results), 1 + length(colnames(GenomicRanges::mcols(wgbs))))
    expect_equal(colnames(results), c("IDs", colnames(GenomicRanges::mcols(wgbs))))
    expect_equal(class(results[[2]]), "numeric")
    expect_equal(class(results[[1]]), "character")
    expect_gt(NROW(results), 0)
    expect_true(all(results$methylation >= 0))
    expect_true(all(results$methylation <= 1))
    expect_true(all(is.element(results$ID, probe_ids$ID)))
    expect_lt(NROW(results), NROW(probe_ids))
    expect_lt(NROW(results), NROW(wgbs))

    expect_error(NROW(BSmeth2Probe(probe_id_locations = probe_ids[0], WGBS_data = wgbs)))
    expect_error(NROW(BSmeth2Probe(probe_id_locations = probe_ids, WGBS_data = wgbs[0])))

    expect_lt(NROW(results), NROW(BSmeth2Probe(probe_id_locations = probe_ids, WGBS_data = wgbs, cutoff = 100)))
    expect_gt(NROW(results), NROW(BSmeth2Probe(probe_id_locations = probe_ids, WGBS_data = wgbs, cutoff = 0)))
    expect_lt(NROW(results), NROW(BSmeth2Probe(probe_id_locations = probe_ids, WGBS_data = wgbs, multipleMapping = TRUE)))
    expect_gt(NROW(results), NROW(BSmeth2Probe(probe_id_locations = probe_ids[1:NROW(probe_ids) / 2], WGBS_data = wgbs)))
    expect_gt(NROW(results), NROW(BSmeth2Probe(probe_id_locations = probe_ids, WGBS_data = wgbs[1:NROW(wgbs) / 2])))

    expect_equal(results, BSmeth2Probe(probe_id_locations = probe_ids_df, WGBS_data = wgbs))
    expect_equal(
        colnames(BSmeth2Probe(probe_id_locations = probe_ids, WGBS_data = methylKit::methRead(list(
            system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
            system.file("extdata", "test2.myCpG.txt", package = "methylKit")
        ),
        sample.id = list("test1", "test2"), assembly = "hg18", treatment = c(1, 1), context = "CpG"
        ))),
        c("IDs", methylKit::getSampleID(methylKit::methRead(list(
            system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
            system.file("extdata", "test2.myCpG.txt", package = "methylKit")
        ),
        sample.id = list("test1", "test2"), assembly = "hg18", treatment = c(1, 1), context = "CpG"
        )))
    )
    expect_equal(
        colnames(BSmeth2Probe(probe_id_locations = probe_ids, WGBS_data = methRead(list(
            system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
            system.file("extdata", "test2.myCpG.txt", package = "methylKit")
        ),
        sample.id = list("test1", "test2"), assembly = "hg18", treatment = c(1, 1), context = "CpG"
        ))),
        colnames(BSmeth2Probe(cutoff = 100, probe_id_locations = probe_ids, WGBS_data = methRead(list(
            system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
            system.file("extdata", "test2.myCpG.txt", package = "methylKit")
        ),
        sample.id = list("test1", "test2"), assembly = "hg18", treatment = c(1, 1), context = "CpG"
        )))
    )

    simulatedData <- methylKit::dataSim(replicates = 4, sites = 200000, treatment = c(1, 1, 0, 0), percentage = 10, effect = 25)
    expect_gte(NROW(BSmeth2Probe(probe_id_locations = probe_ids_df, WGBS_data = simulatedData)), 0)
    expect_equal(
        NCOL(BSmeth2Probe(probe_id_locations = probe_ids_df, WGBS_data = simulatedData)),
        length(methylKit::getSampleID(simulatedData)) + 1
    )


    expect_error(BSmeth2Probe(probe_id_locations = probe_ids, WGBS_data = wgbs, cutoff = -1))
    expect_error(BSmeth2Probe(probe_id_locations = probe_ids, WGBS_data = probe_ids))
    expect_error(BSmeth2Probe(probe_id_locations = wgbs, WGBS_data = wgbs))
    expect_error(BSmeth2Probe(probe_id_locations = NULL, WGBS_data = wgbs))
    expect_error(BSmeth2Probe(probe_id_locations = probe_ids, WGBS_data = NULL))
})
