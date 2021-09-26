
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
User who wish to expand their own reference atlas can use `findSignatures` function. Below is an example dataframe to illustrate the `atlas` format. `atlas` is the signature matrix to be extended and `samples` the new data to be added to the signature matrix. `atlas` and `samples` are compliant with the function requirements. After providing appropriate `atlas` format, users can create `samples` using `simulateCellMix` function.

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

