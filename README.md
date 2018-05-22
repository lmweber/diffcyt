# diffcyt

[![Build Status](https://travis-ci.org/lmweber/diffcyt.svg?branch=master)](https://travis-ci.org/lmweber/diffcyt)


## Summary

`diffcyt`: R package for differential discovery in high-dimensional cytometry via high-resolution clustering

The `diffcyt` package implements statistical methods for differential discovery analyses in high-dimensional cytometry data (including flow cytometry, mass cytometry or CyTOF, and oligonucleotide-tagged cytometry), based on (i) high-resolution clustering and (ii) moderated tests adapted from transcriptomics.


## Details

For details on the statistical methodology and comparisons with existing approaches, see the accompanying paper (available soon).


## Tutorial and examples

For a tutorial and examples of usage, see the package vignette.


## Availability and installation

The `diffcyt` package is available from [Bioconductor](http://bioconductor.org/packages/diffcyt). It can be installed using the Bioconductor installer (`biocLite`):

```{r}
# Download the Bioconductor installer
source("https://bioconductor.org/biocLite.R")

# Install 'diffcyt' package
biocLite("diffcyt")
```


The development version of the `diffcyt` package is available from the `devel` version of Bioconductor or from GitHub. See the Bioconductor help pages for how to set up the `devel` version of Bioconductor.

To install the development version from GitHub, the `devtools` package can be used. Note that dependency packages will need to be installed separately from Bioconductor and CRAN.

```{r}
# Install and load 'devtools' package
install.packages("devtools")
library(devtools)

# Install development version of 'diffcyt' from GitHub
install_github("lmweber/diffcyt")
```

