# diffcyt

[![R build status](https://github.com/lmweber/diffcyt/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/lmweber/diffcyt/actions)


## Introduction

`diffcyt`: R package for differential discovery in high-dimensional cytometry via high-resolution clustering

The `diffcyt` package implements statistical methods for differential discovery analyses in high-dimensional cytometry data (including flow cytometry, mass cytometry or CyTOF, and oligonucleotide-tagged cytometry), based on a combination of high-resolution clustering and empirical Bayes moderated tests adapted from transcriptomics.

<p> <img src="vignettes/diffcyt.png" width="130"/> </p>


## Details and citation

For details on the statistical methodology and comparisons with existing approaches, see our paper:

- [Weber et al. (2019), *diffcyt: Differential discovery in high-dimensional cytometry via high-resolution clustering*, Communications Biology, 2, 183](https://www.nature.com/articles/s42003-019-0415-5)


## Tutorial and examples

For a tutorial and examples of usage, see the Bioconductor [package vignette](http://bioconductor.org/packages/release/bioc/vignettes/diffcyt/inst/doc/diffcyt_workflow.html) (link also available via the main Bioconductor page for the [diffcyt package](http://bioconductor.org/packages/diffcyt)).


## Installation

The `diffcyt` package is available from [Bioconductor](http://bioconductor.org/packages/diffcyt), and can be installed as follows:

```{r}
# Install Bioconductor installer from CRAN
install.packages("BiocManager")

# Install 'diffcyt' package from Bioconductor
BiocManager::install("diffcyt")
```

To run the examples in the package vignette and generate additional visualizations, the `HDCytoData` and `CATALYST` packages from Bioconductor are also required.

```{r}
BiocManager::install(c("HDCytoData", "CATALYST"))
```
