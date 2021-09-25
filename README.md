
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

**deconvR** is a collection of functions surrounding the deconvolution of bulk
sample(s) with the use of a reference atlas of signature profiles and a
user-selected model. The users can use create/extend their reference atlas using
`findSignatures` function or, they may choose to use the reference atlas provided within the package.

## Installation

The deconvR package can be installed from Bioconductor with:

``` {r }
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("deconvR")
```
You can also install deconvR directly with github:

``` r

devtools::install_github("BIMSBbioinfo/deconvR")
```
