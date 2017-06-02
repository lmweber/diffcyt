# diffcyt

R package: Statistical methods for differential discovery in high-dimensional flow cytometry and mass cytometry (CyTOF) data.

Under development.



## How to install

Bioconductor dependencies need to be installed manually (since `install_github` can only install CRAN dependencies automatically):

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite(c("flowCore", "FlowSOM", "limma", "IHW", "SummarizedExperiment", "BiocParallel"))
```

After the Bioconductor dependencies are installed, the package can be installed from GitHub. Note that an authentication token is required since this is a private repository.

```{r}
library(devtools)
install_github("lmweber/diffcyt", auth_token = "...")
```



## Types of differential discovery

We consider two main types of differential discovery / analysis:

- **Differential abundance** of cell populations. Cell populations are defined using high-resolution clustering.

- **Differential expression of functional markers** within cell populations. Cell populations are defined using high-resolution clustering on a subset of protein markers (e.g. surface markers in immunology). Additional functional markers are then analyzed for differential expression within the clusters. The two sets of markers (clustering and functional) must be specified by the user.



## Statistical approaches

Several different statistical approaches are under development:

- `diffcyt-DA`: differential abundance of clusters; based on `limma-voom` methods.

- `diffcyt-med`: differential expression of functional markers within clusters; simple approach based on medians for each functional marker expression profile

- `diffcyt-FDA`: differential expression of functional markers within clusters; using methods from functional data analysis (FDA) to test for differences between functional marker expression profiles represented by empirical cumulative distribution functions (ECDFs) in each condition

- `diffcyt-KS`: differential expression of functional markers within clusters; using Kolmogorov-Smirnov tests to test for differences between ECDFs in each condition

- `diffcyt-LM`: differential expression of functional markers within clusters; using linear models or mixed models to test for differences between ECDFs in each condition


