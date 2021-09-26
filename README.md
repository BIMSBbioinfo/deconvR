
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

The **deconvR** is an [R](http://en.wikipedia.org/wiki/R_%28programming_language%29) package for analyzing deconvolution of the bulk sample(s) using an atlas of reference signature profiles and a user-selected model. Users can create or expand their own reference atlases using the `findSignatures` function, or they can choose to use the reference atlas provided in the package. A more detailed explanation on to use deconvR can be found in `How to Use deconvR` section.

## Installation

The deconvR package can be installed from Bioconductor with:

``` {r }
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("deconvR")
```
You can also install deconvR directly with github:

``` r
remotes::install_github("BIMSBbioinfo/deconvR")
```

## How to Use deconvR


