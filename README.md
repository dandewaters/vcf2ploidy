# vcf2ploidy 2.0

## Introduction

This package is an update to my original [VCF File Converter tool](https://github.com/dandewaters/VCF-File-Converter), rewritten for R. This R package acts as a companion to the [gbs2ploidy package](https://cran.r-project.org/package=gbs2ploidy). The functions in the gbs2ploidy package require the input data to be in Heterozygous Allele Depth (HAD) format, however GBS2Ploidy does not include a function for converting from Variant-Calling Format (VCF). The VCF2HAD() function in this package reads in VCF files and returns a data frame in HAD format that can be input to the GBS2Ploidy functions. Alternatively, you can use the vcf2ploidy() function in this package to achieve this in one step.

This package also contains a function for converting VCF files to Colony format.

## Changes with version 2.0!

* Function names are now all lowercase for ease of use and to match gbs2ploidy's style
* vcf2ploidy now has an opional GUI using shiny!

  Simply run the command `vcf2ploidy_app()` and the interface will launch in RStudio's "Viewer" panel:

![](inst/vcf2ploidy_app.png)


## Installation
 
To install, make sure you have the [devtools package](https://cran.r-project.org/package=devtools) installed and loaded. Then run the following commands:

```{r installation, eval=FALSE}
library(devtools)
install_github("dandewaters/vcf2ploidy")
```

## Vignettes

Read the intro vignette by running install_github with vignettes = TRUE and running the following commands:

```{r vignettes, eval=FALSE}
library(devtools)
install_github("dandewaters/vcf2ploidy", build_vignette=TRUE)
vignette("Introduction", package="vcf2ploidy")
```

## How to Cite Me
CSE Bibliography Format:

DeWaters, D. 2020. vcf2ploidy. Waldorf (MD): GitHub; [accessed Year Month Day]. https://github.com/dandewaters/vcf2ploidy.

## References

Gompert  Z.  and  Mock  K.  (XXXX)  Detection  of  individual  ploidy  levels  with  genotyping-by-sequencing (GBS) analysis. Molecular Ecology Resources, submitted.
Jones OR, Wang J. COLONY: a program for parentage and sibship inference from multilocus genotype data. Mol Ecol Resour. 2010 May;10(3):551-5. doi: 10.1111/j.1755-0998.2009.02787.x. Epub 2009 Oct 21. PMID: 21565056.
