
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GSAMU

<!-- badges: start -->

[![Project
Status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active/)
[![Package
version](https://img.shields.io/badge/GitHub-1.0.0-orange.svg)](https://github.com/seungjae2525/GSAMU/)
[![minimal R
version](https://img.shields.io/badge/R-v4.0.0+-blue.svg)](https://cran.r-project.org/)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![R-CMD-check](https://github.com/seungjae2525/GSAMU/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/seungjae2525/GSAMU/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Description

This is the source code for the `GSAMU` package in R. `GSAMU` is a
package aimed at providing a novel sensitivity model to investigate the
effect of correlated multiple exposures on the non-Gaussian and
time-to-event outcome variables. Given a user-specified sensitivity
parameters, the sensitivity range is calculated. See reference for
details.

## Reference

Lee S, Jeong B, Lee D, Lee W (2024): Sensitivity analysis for effects of
multiple exposures in the presence of unmeasured confounding:
non-Gaussian and time-to-event outcomes. submitted.

## Installation

### Current GitHub release:

Installation using R package `remotes`:

``` r
install.packages("remotes") # if devtools not already installed
remotes::install_github("seungjae2525/GSAMU")
library(GSAMU)
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(GSAMU)
### basic example code
## For count outcome, 
count.re <- GSAMU(data=dataset, 
                  outcome="Y", outcome.type="count", link="log",
                  confounder=c("L1", "L2", "L3"),
                  exposure=c("X1", "X2", "X3", "X4"),
                  delta.range=c(0.11, 0.44), delta.diff=0.11, k=k, p=p, bound=bound,
                  bootsCI=FALSE, B=1000, seed=231111, verbose=TRUE,
                  report.result=TRUE, decimal.p=3)

## For binary outcome with logit link, 
binary.re1 <- GSAMU(data=dataset, 
                    outcome="Y", outcome.type="count", link="log",
                    confounder=c("L1", "L2", "L3"),
                    exposure=c("X1", "X2", "X3", "X4"),
                    delta.range=c(0.11, 0.44), delta.diff=0.11, k=k, p=p, bound=bound,
                    bootsCI=FALSE, B=1000, seed=231111, verbose=TRUE,
                    report.result=TRUE, decimal.p=3)

## For binary outcome with logit link, 
binary.re2 <-GSAMU(data=dataset, 
                   outcome="Y", outcome.type="count", link="log",
                   confounder=c("L1", "L2", "L3"),
                   exposure=c("X1", "X2", "X3", "X4"),
                   delta.range=c(0.11, 0.44), delta.diff=0.11, k=k, p=p, bound=bound,
                   bootsCI=FALSE, B=1000, seed=231111, verbose=TRUE,
                   report.result=TRUE, decimal.p=3)

## For time-to-event outcome, 
cox.re <- GSAMU(data=dataset, 
                outcome=c("time", "status"), outcome.type="timetoevent", link=NULL,
                confounder=c("L1", "L2", "L3"),
                exposure=c("X1", "X2", "X3", "X4"),
                delta.range=c(0.11, 0.44), delta.diff=0.11, k=k, p=p, bound=bound,
                bootsCI=FALSE, B=1000, seed=231111, verbose=TRUE,
                report.result=TRUE, decimal.p=3)
```

You can also resulted plots, for example:

``` r
autoplot(object=count.re, point.size=2.75, width.SI=1.55, width.CI=0.6,
         axis.title.x.size=15, axis.text.size=16, legend.text.size=15,
         myxlim=c(-0.25, 2))

autoplot(object=binary.re1, point.size=2.75, width.SI=1.55, width.CI=0.6,
         axis.title.x.size=15, axis.text.size=16, legend.text.size=15,
         myxlim=c(-0.25, 2))

autoplot(object=binary.re2, point.size=2.75, width.SI=1.55, width.CI=0.6,
         axis.title.x.size=15, axis.text.size=16, legend.text.size=15,
         myxlim=c(-0.25, 2))

autoplot(object=cox.re, point.size=2.75, width.SI=1.55, width.CI=0.6,
         axis.title.x.size=15, axis.text.size=16, legend.text.size=15,
         myxlim=c(-0.25, 2))
```

## Bug Reports:

You can also report bugs on GitHub under
[Issues](https://github.com/seungjae2525/GSAMU/issues/).
