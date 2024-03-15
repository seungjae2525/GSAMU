# GSAMU

<!-- badges: start -->
[![Project Status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active/)
[![Package version](https://img.shields.io/badge/GitHub-1.0.0-orange.svg)](https://github.com/seungjae2525/GSAMU/)
[![minimal R version](https://img.shields.io/badge/R-v4.0.0+-blue.svg)](https://cran.r-project.org/)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![R-CMD-check](https://github.com/seungjae2525/GSAMU/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/seungjae2525/GSAMU/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# Sensitivity analysis for effects of multiple exposures in the presence of unmeasured confounding: non-Gaussian and time-to-event outcomes

## Description
This is the source code for the `GSAMU` package in R. 
`GSAMU` is a package aimed at providing a novel sensitivity model to investigate the effect of correlated multiple exposures on the non-Gaussian and time-to-event outcome variables.
Given a user-specified sensitivity parameters, the sensitivity range is calculated. See reference for details.

### Reference
Lee S, Jeong B, Lee D, Lee W (2024): Sensitivity analysis for effects of multiple exposures in the presence of unmeasured confounding: non-Gaussian and time-to-event outcomes. submitted.


## Installation
### Current GitHub release:
Installation using R package `remotes`:

```r
install.packages("remotes") # if devtools not already installed
remotes::install_github("seungjae2525/GSAMU")
library(GSAMU)
```

### Bug Reports:
You can also report bugs on GitHub under [Issues](https://github.com/seungjae2525/GSAMU/issues/).
