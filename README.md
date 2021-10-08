
<!-- README.md is generated from README.Rmd. Please edit that file -->

# deconvR : Simulation and Deconvolution of Omic Profiles <img src="deconvR_logo.png" align="right"  alt="logo" width="250" style = "border: none; float: left ;">

[![R-CMD-check-bioc](https://github.com/BIMSBbioinfo/deconvR/actions/workflows/check-bioc.yml/badge.svg)](https://github.com/BIMSBbioinfo/deconvR/actions/workflows/check-bioc.yml)
[![R-CMD-check](https://github.com/BIMSBbioinfo/deconvR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/BIMSBbioinfo/deconvR/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/BIMSBbioinfo/deconvR/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/BIMSBbioinfo/deconvR/actions/workflows/test-coverage.yaml)
[![codecov](https://codecov.io/gh/BIMSBbioinfo/deconvR/branch/master/graph/badge.svg?token=F86XU6BI9S)](https://codecov.io/gh/BIMSBbioinfo/deconvR)
[![BioCstatus](http://www.bioconductor.org/shields/build/release/bioc/deconvR.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/deconvR)
[![](https://img.shields.io/badge/download-NA/total-blue.svg)](https://bioconductor.org/packages/stats/bioc/deconvR)

<!-- badges: start -->
<!-- badges: end -->

The **deconvR** package designed for analyzing deconvolution of the bulk
sample(s) using an atlas of reference signature profiles and a
user-selected model (non-negative least squares regression, support
vector regression, quadratic programming, or robust linear regression).
Users can upload or expand their own reference atlases using the
`findSignatures` function, or they can choose to use the reference atlas
provided in the package. `simulateCellMix` function included to simulate
a bulk signature profile of a given size. Additionnaly, `BSmeth2Probe`
function can be used to map methylation data to probe IDs.

## Installation

The deconvR package can be installed from Bioconductor with:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("deconvR")
```

You can also install the development version of the **deconvR** directly
from GitHub:

``` r
remotes::install_github("BIMSBbioinfo/deconvR")
```

## How to Use deconvR

User who wish to expand their own reference atlas can use
`findSignatures` function. `atlas` is the signature matrix to be
extended and `samples` the new data to be added to the signature matrix.
`atlas` and `samples` are compliant with the function requirements.
After providing appropriate `atlas` format, users can create `samples`
using `simulateCellMix` function.
