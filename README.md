# VCF2Ploidy

## Introduction

This package is an update to my original [VCF File Converter tool](https://github.com/dandewaters/VCF-File-Converter), rewritten for R. This R package acts as a companion to the [gbs2ploidy package](https://cran.r-project.org/package=gbs2ploidy). The functions in the GBS2Ploidy package require the input data to be in Heterozygous Allele Depth (HAD) format, however GBS2Ploidy does not include a function for converting from Variant-Calling Format (VCF). The VCF2HAD() function in this package reads in VCF files and returns a data frame in HAD format that can be input to the GBS2Ploidy functions. Alternatively, you can use the VCF2Ploidy() function in this package to achieve this in one step.

## Installation
 
To install, make sure you have the [devtools package](https://cran.r-project.org/package=devtools) installed and loaded. Then run the following commands:

```{r installation, eval=FALSE}
library(devtools)
install_github("dandewaters/VCF2Ploidy")
library(FARSFunctions)
```

## Vignettes

Read the intro vignette by running install_github with vignettes = TRUE and running the following commands:

```{r vignettes, eval=FALSE}
library(devtools)
install_github("dandewaters/VCF2Ploidy", build_vignette=TRUE)
vignette("Introduction", package="VCF2Ploidy")
```

<!-- badges: start -->
[![Travis build status](https://travis-ci.com/dandewaters/VCF2Ploidy.svg?branch=master)](https://travis-ci.com/dandewaters/VCF2Ploidy)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/dandewaters/VCF2Ploidy?branch=master&svg=true)](https://ci.appveyor.com/project/dandewaters/VCF2Ploidy)

<!-- badges: end -->
