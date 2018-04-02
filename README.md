# diffcyt

[![Build Status](https://travis-ci.org/lmweber/diffcyt.svg?branch=master)](https://travis-ci.org/lmweber/diffcyt)


## Summary

`diffcyt`: R package for differential discovery in high-dimensional cytometry via high-resolution clustering

The `diffcyt` package implements statistical methods for differential discovery analyses in high-dimensional cytometry data (including flow cytometry, mass cytometry or CyTOF, and DNA-tagged cytometry), based on high-resolution clustering and moderated tests.


## Details

For details on the statistical methodology and comparisons with existing approaches, see the accompanying paper (available soon).


## Tutorial and examples

For a tutorial and examples of usage, see the package vignette.


## Availability and installation

The stable release version of the `diffyt` package will be made available from [Bioconductor](http://bioconductor.org/).

The development version is available from the `bioc-devel` version of Bioconductor, or from GitHub.

To install from GitHub, the Bioconductor installer (`biocLite`) can be used. This will also install all required dependencies from CRAN and Bioconductor. First, ensure that the Bioconductor installer and `devtools` package are installed:

```{r}
source("https://bioconductor.org/biocLite.R")
install.packages("devtools")
```

Then, install the `diffcyt` package from GitHub using `biocLite`, with the option `dependencies = TRUE`.

```{r}
biocLite("lmweber/diffcyt", dependencies = TRUE)
```


