
<!-- README.md is generated from README.Rmd. Please edit that file -->

# deconvR
[![Build Status](https://app.travis-ci.com/BIMSBbioinfo/deconvR.svg?branch=main)](https://app.travis-ci.com/github/BIMSBbioinfo/deconvR)
<!-- badges: start -->
<!-- badges: end -->

The package contains the function “deconvolute” which takes in a
reference atlas (defaults to ref\_atlas included in package), a bulk
sample(s) signature profile, and the desired model to be used
(non-negative least squares regression, support vector regression,
quadratic programming, or robust linear regression). The desired model
will be used in conjunction with the reference atlas in order to predict
the cell origin profiles of the bulk sample(s). The partial R-squared
values of the model are also printed upon completion. Also included is
the function simulateCellMix, which can be used to simulate a bulk
signature profile of a given size, using a given reference atlas (again
defaults to ref\_atlas included in package). The outputted table of bulk
samples is compatible to be used as the bulk input for deconvolute. Also
included is BSmeth2Probe, which can be used to map WGBS methylation data
to probe IDs. The output dataframe of BSmeth2Probe can be used for the
bulk sample methylation profile of deconvolute. Also included is
findSignatures, which can create or extend a reference atlas containing
signatures of cell types.

## Installation

You can install deconvR with:

``` r
devtools::install_github("BIMSBbioinfo/deconvR")
```
