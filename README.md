
<!-- README.md is generated from README.Rmd. Please edit that file -->
<a name="deconvR_logo"/>
<div align="center">
<img src="https://github.com/BIMSBbioinfo/deconvR/blob/main/inst/deconvR_logo.png" alt="deconvR_logo" width="650"/ ></img>
</a>
</div>

# deconvR : Simulation and Deconvolution of Cellular Signatures
[![R-CMD-check-bioc](https://github.com/BIMSBbioinfo/deconvR/actions/workflows/check-bioc.yml/badge.svg)](https://github.com/BIMSBbioinfo/deconvR/actions/workflows/check-bioc.yml)   [![R-CMD-check](https://github.com/BIMSBbioinfo/deconvR/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/BIMSBbioinfo/deconvR/actions/workflows/check-standard.yaml) [![test-coverage](https://github.com/BIMSBbioinfo/deconvR/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/BIMSBbioinfo/deconvR/actions/workflows/test-coverage.yaml)  [![codecov](https://codecov.io/gh/BIMSBbioinfo/deconvR/branch/main/graph/badge.svg)](https://github.com/BIMSBbioinfo/deconvR/actions)   [![BioC status](http://www.bioconductor.org/shields/build/release/bioc/deconvR.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/deconvR)


<!-- badges: start -->
<!-- badges: end -->

The **deconvR** is an [R](http://en.wikipedia.org/wiki/R_%28programming_language%29) package for analyzing deconvolution of the bulk sample(s) using an atlas of reference signature profiles and a user-selected model. Users can upload or expand their own reference atlases using the `findSignatures` function, or they can choose to use the reference atlas provided in the package. A more detailed explanation on to use deconvR can be found in `How to Use deconvR` section.

## Installation

The deconvR package can be installed from Bioconductor with:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("deconvR")
```
You can also install deconvR directly with github:

``` r
remotes::install_github("BIMSBbioinfo/deconvR")
```

## How to Use deconvR
#### Expanding Reference Atlas
User who wish to expand their own reference atlas can use `findSignatures` function. Below is an example dataframe to illustrate the `atlas` format. `atlas` is the signature matrix to be extended and `samples` the new data to be added to the signature matrix. `atlas` and `samples` are compliant with the function requirements.

``` r 
atlas <- readRDS(system.file("reference_atlas_nodup.RDS", package = "deconvR")) #Read the reference atlas provided within the package
```
Such `atlas` format must provided in order to use `findSignatures` function.
``` r 
        IDs Monocytes_EPIC B.cells_EPIC CD4T.cells_EPIC NK.cells_EPIC CD8T.cells_EPIC
1 cg08169020         0.8866       0.2615          0.0149        0.0777          0.0164
2 cg25913761         0.8363       0.2210          0.2816        0.4705          0.3961
3 cg26955540         0.7658       0.0222          0.1492        0.4005          0.3474
4 cg25170017         0.8861       0.5116          0.1021        0.4363          0.0875
5 cg12827637         0.5212       0.3614          0.0227        0.2120          0.0225
6 cg19442545         0.2013       0.1137          0.0608        0.0410          0.0668
```
`samples` can be created using `simulateCellMix` function.

``` r 
samples <- simulateCellMix(3)[[1]]
```

The format of the sample is as follows:

``` r 
         IDs  Sample 1   Sample 2 Sample 3
1 cg08169020 0.1246110 0.05314110   0.0168
2 cg25913761 0.2923030 0.28962997   0.3104
3 cg26955540 0.1717311 0.15979721   0.0978
4 cg25170017 0.2500497 0.14630151   0.2832
5 cg12827637 0.1565242 0.08521324   0.1368
6 cg19442545 0.0589429 0.05893091   0.0222
```

`sampleMeta` should include all sample names in `samples`, and specify the origins they should be mapped to when added to `atlas`.

``` r 
sampleMeta <- data.table("Experiment_accession" = colnames(samples)[-1], "Biosample_term_name" = "new cell type")
```

Example of `sampleMeta` as follows:
```r
   Experiment_accession Biosample_term_name
1:             Sample 1       new cell type
2:             Sample 2       new cell type
3:             Sample 3       new cell type
```

Then we can use `findSignatures` to extend the matrix.
``` r
extended_matrix <- findSignatures(samples = samples, 
                                 sampleMeta = sampleMeta, 
                                 atlas = atlas)
```

Output format of the extended matrix obtained as follows:
```r
          IDs new_cell_type Monocytes_EPIC B.cells_EPIC CD4T.cells_EPIC NK.cells_EPIC
1: cg08169020    0.06485071         0.8866       0.2615          0.0149        0.0777
2: cg25913761    0.29744433         0.8363       0.2210          0.2816        0.4705
3: cg26955540    0.14310944         0.7658       0.0222          0.1492        0.4005
4: cg25170017    0.22651708         0.8861       0.5116          0.1021        0.4363
5: cg12827637    0.12617913         0.5212       0.3614          0.0227        0.2120
6: cg19442545    0.04669127         0.2013       0.1137          0.0608        0.0410
```
