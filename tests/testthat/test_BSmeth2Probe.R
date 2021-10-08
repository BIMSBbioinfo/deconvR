library(deconvR)

test_that("BSmeth2Probe", {
    WGBS_GRanges <- readRDS(system.file("extdata", "WGBS_GRanges.RDS",
        package = "deconvR"
    ))
    data("IlluminaInfiniumMethylationEpicv1B5ManifestProbes")
    IlluminaInfiniumMethylationEpicv1B5ManifestProbes_df <-
        GenomicRanges::as.data.frame(IlluminaInfiniumMethylationEpicv1B5ManifestProbes)

    results <- BSmeth2Probe(
        probe_id_locations = IlluminaInfiniumMethylationEpicv1B5ManifestProbes,
        WGBS_data = WGBS_GRanges
    )

    expect_equal(class(results), "data.frame")
    expect_equal(NCOL(results), 1 +
        length(colnames(GenomicRanges::mcols(WGBS_GRanges))))
    expect_equal(colnames(results), c(
        "IDs",
        colnames(GenomicRanges::mcols(WGBS_GRanges))
    ))
    expect_equal(class(results[[2]]), "numeric")
    expect_equal(class(results[[1]]), "character")
    expect_gt(NROW(results), 0)
    expect_true(all(results$methylation >= 0))
    expect_true(all(results$methylation <= 1))
    expect_true(all(is.element(
        results$ID,
        IlluminaInfiniumMethylationEpicv1B5ManifestProbes$ID
    )))
    expect_lt(
        NROW(results),
        NROW(IlluminaInfiniumMethylationEpicv1B5ManifestProbes)
    )
    expect_lt(NROW(results), NROW(WGBS_GRanges))

    expect_error(NROW(BSmeth2Probe(
        probe_id_locations =
            IlluminaInfiniumMethylationEpicv1B5ManifestProbes[0],
        WGBS_data = WGBS_GRanges
    )))
    expect_error(NROW(BSmeth2Probe(
        probe_id_locations = IlluminaInfiniumMethylationEpicv1B5ManifestProbes,
        WGBS_data = WGBS_GRanges[0]
    )))

    expect_lt(NROW(results), NROW(BSmeth2Probe(
        probe_id_locations = IlluminaInfiniumMethylationEpicv1B5ManifestProbes,
        WGBS_data = WGBS_GRanges,
        cutoff = 100
    )))
    expect_gt(NROW(results), NROW(BSmeth2Probe(
        probe_id_locations = IlluminaInfiniumMethylationEpicv1B5ManifestProbes,
        WGBS_data = WGBS_GRanges,
        cutoff = 0
    )))
    expect_lt(NROW(results), NROW(BSmeth2Probe(
        probe_id_locations = IlluminaInfiniumMethylationEpicv1B5ManifestProbes,
        WGBS_data = WGBS_GRanges,
        multipleMapping = TRUE
    )))
    expect_gt(NROW(results), NROW(BSmeth2Probe(
        probe_id_locations =
            IlluminaInfiniumMethylationEpicv1B5ManifestProbes
            [1:NROW(IlluminaInfiniumMethylationEpicv1B5ManifestProbes) / 2],
        WGBS_data = WGBS_GRanges
    )))
    expect_gt(NROW(results), NROW(BSmeth2Probe(
        probe_id_locations = IlluminaInfiniumMethylationEpicv1B5ManifestProbes,
        WGBS_data =
            WGBS_GRanges[1:NROW(WGBS_GRanges) / 2]
    )))

    expect_equal(results, BSmeth2Probe(
        probe_id_locations =
            IlluminaInfiniumMethylationEpicv1B5ManifestProbes_df,
        WGBS_data = WGBS_GRanges
    ))
    expect_equal(
        colnames(BSmeth2Probe(
            probe_id_locations =
                IlluminaInfiniumMethylationEpicv1B5ManifestProbes,
            WGBS_data = methylKit::methRead(list(
                system.file("extdata", "test1.myCpG.txt",
                    package = "methylKit"
                ),
                system.file("extdata", "test2.myCpG.txt",
                    package = "methylKit"
                )
            ),
            sample.id = list("test1", "test2"), assembly = "hg18",
            treatment = c(1, 1), context = "CpG"
            )
        )),
        c("IDs", methylKit::getSampleID(methylKit::methRead(list(
            system.file("extdata", "test1.myCpG.txt",
                package = "methylKit"
            ),
            system.file("extdata", "test2.myCpG.txt",
                package = "methylKit"
            )
        ),
        sample.id = list("test1", "test2"), assembly = "hg18",
        treatment = c(1, 1), context = "CpG"
        )))
    )
    expect_equal(
        colnames(BSmeth2Probe(
            probe_id_locations =
                IlluminaInfiniumMethylationEpicv1B5ManifestProbes,
            WGBS_data = methylKit::methRead(list(
                system.file("extdata", "test1.myCpG.txt",
                    package = "methylKit"
                ),
                system.file("extdata", "test2.myCpG.txt",
                    package = "methylKit"
                )
            ),
            sample.id = list("test1", "test2"), assembly = "hg18",
            treatment = c(1, 1), context = "CpG"
            )
        )),
        colnames(BSmeth2Probe(
            cutoff = 100, probe_id_locations =
                IlluminaInfiniumMethylationEpicv1B5ManifestProbes,
            WGBS_data = methylKit::methRead(list(
                system.file("extdata", "test1.myCpG.txt",
                    package = "methylKit"
                ),
                system.file("extdata", "test2.myCpG.txt",
                    package = "methylKit"
                )
            ),
            sample.id = list("test1", "test2"), assembly = "hg18",
            treatment = c(1, 1), context = "CpG"
            )
        ))
    )

    simulatedData <- methylKit::dataSim(
        replicates = 4, sites = 200000,
        treatment = c(1, 1, 0, 0),
        percentage = 10, effect = 25
    )
    expect_gte(NROW(BSmeth2Probe(
        probe_id_locations =
            IlluminaInfiniumMethylationEpicv1B5ManifestProbes_df,
        WGBS_data = simulatedData
    )), 0)
    expect_equal(
        NCOL(BSmeth2Probe(
            probe_id_locations =
                IlluminaInfiniumMethylationEpicv1B5ManifestProbes_df,
            WGBS_data = simulatedData
        )),
        length(methylKit::getSampleID(simulatedData)) + 1
    )


    expect_error(BSmeth2Probe(
        probe_id_locations =
            IlluminaInfiniumMethylationEpicv1B5ManifestProbes,
        WGBS_data = WGBS_GRanges, cutoff = -1
    ))
    expect_error(BSmeth2Probe(
        probe_id_locations =
            IlluminaInfiniumMethylationEpicv1B5ManifestProbes,
        WGBS_data = IlluminaInfiniumMethylationEpicv1B5ManifestProbes
    ))
    expect_error(BSmeth2Probe(
        probe_id_locations = WGBS_GRanges,
        WGBS_data = WGBS_GRanges
    ))
    expect_error(BSmeth2Probe(
        probe_id_locations = NULL,
        WGBS_data = WGBS_GRanges
    ))
    expect_error(BSmeth2Probe(
        probe_id_locations =
            IlluminaInfiniumMethylationEpicv1B5ManifestProbes,
        WGBS_data = NULL
    ))
})
