# diffcyt

[![Build Status](https://travis-ci.org/lmweber/diffcyt.svg?branch=master)](https://travis-ci.org/lmweber/diffcyt)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/lmweber/diffcyt?branch=master&svg=true)](https://ci.appveyor.com/project/lmweber/diffcyt)

R package: Differential discovery in high-dimensional cytometry via high-resolution clustering


## Details

Statistical methods for differential discovery in high-dimensional cytometry (including flow cytometry, mass cytometry or CyTOF, and DNA-tagged cytometry) using high-resolution clustering and moderated tests.

For details on the methodology and comparisons with existing approaches, see the accompanying paper (available soon).


## Tutorial and examples

For a tutorial and examples of usage, see the vignette.


## How to install

The package will be submitted to Bioconductor (soon). It can also be installed from GitHub.

The Bioconductor installer (`biocLite`) can be used to install from GitHub. This will install all dependencies from both CRAN and Bioconductor.

First, ensure that the Bioconductor installer and `devtools` package are installed:

```{r}
source("https://bioconductor.org/biocLite.R")
install.packages("devtools")
```

Then, install the `diffcyt` package from GitHub using `biocLite`, with the option `dependencies = TRUE`.

```{r}
biocLite("lmweber/diffcyt", dependencies = TRUE)
```


