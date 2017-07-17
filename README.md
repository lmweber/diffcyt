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

We consider two main types of differential discovery:

- **Differential abundance** of cell populations. Cell populations are defined using high-resolution clustering.

- **Differential expression of functional markers** within cell populations. Cell populations are defined using high-resolution clustering on a subset of protein markers (e.g. surface markers in immunology). Additional functional markers are then analyzed for differential expression within the clusters. The two sets of markers (clustering and functional) must be specified by the user.



## Statistical approaches

Several different statistical approaches are under development. More details later.


