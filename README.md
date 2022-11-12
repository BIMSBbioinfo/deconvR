
<!-- README.md is generated from README.Rmd. Please edit that file -->

# deconvR : Simulation and Deconvolution of Omic Profiles

[![R-CMD-check](https://github.com/BIMSBbioinfo/deconvR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/BIMSBbioinfo/deconvR/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/BIMSBbioinfo/deconvR/branch/master/graph/badge.svg?token=F86XU6BI9S)](https://codecov.io/gh/BIMSBbioinfo/deconvR)
[![](https://img.shields.io/badge/release%20version-1.4.1-green.svg)](https://www.bioconductor.org/packages/deconvR)

<!-- badges: start -->
<!-- badges: end -->

<img src="deconvR_logo.png" align="left" alt="logo" width="300" style = "border: none; float: center ;">

The **deconvR** package provides a collection of functions designed for
analyzing deconvolution of the bulk sample(s) using an atlas of
reference signature profiles and a user-selected model (non-negative
least squares,quadratic programming, support vector regression, or
robust linear regression). Users can directly use their reference atlas
or, create an expended version of their reference atlas using
`findSignatures`. Additionnaly, they can also use the reference atlas
provided within the package, which contains cell-type specific
methylation values. Another option is to simulate bulk signatures of
bulk samples using `simulateCellMix`. And finally, we included
`BSmeth2Probe` function along with the `Illumina Methylation EPIC B5 Manifest`
probe IDs, to be used to map methylation data to respective probe IDs.

## Installation

The deconvR package can be installed from Bioconductor:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("deconvR")
```

You can also install the development version of **deconvR** directly
from Github:

``` r
remotes::install_github("BIMSBbioinfo/deconvR")
```

## How to Use deconvR

User who wish to expand their own reference atlas can use
`findSignatures` function. `atlas` is the signature matrix to be
extended and `samples` the new data to be added to the signature matrix.
`atlas` and `samples` are compliant with the function requirements.
After providing appropriate `atlas` format, users can create `samples`
using `simulateCellMix` function. You can get more information about
**deconvR** from [here.](http://bioinformatics.mdc-berlin.de/deconvR/)
