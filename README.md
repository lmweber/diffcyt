# diffcyt

R package: Statistical methods for differential discovery in high-dimensional cytometry data.

Under development.



## Details and examples

For details on the statistical methodology and comparisons with existing approaches, see the paper (in progress).

For a tutorial and examples of usage, see the Bioconductor vignette (in progress).



## How to install

Since there are dependencies from both CRAN and Bioconductor, the Bioconductor installer (`biocLite`) should be used. This will install all dependencies automatically.

Ensure the Bioconductor installer and `devtools` are installed:

```{r}
source("https://bioconductor.org/biocLite.R")
install.packages("devtools")
```

Install the `diffcyt` package from GitHub using `biocLite`, with the option `dependencies = TRUE`. (Currently, an authentication token is also required with `auth_token = ...`, since this is a private repository.)

```{r}
biocLite("lmweber/diffcyt", dependencies = TRUE, auth_token = "...")
```


