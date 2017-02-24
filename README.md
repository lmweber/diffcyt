# diffcyt

R package: Statistical methods for differential discovery in high-dimensional flow cytometry and mass cytometry (CyTOF) data.

Work in progress.



## Types of differential discovery

Methods are under development for two main types of differential discovery analysis:

- *Differential abundance* of cell populations, using [FlowSOM](http://bioconductor.org/packages/release/bioc/html/FlowSOM.html) for automated clustering, and empirical Bayes methods to improve power by sharing information on variability across clusters.

- *Differential expression of functional markers* within automatically defined clusters. For example, in immunology, cell populations (clusters) may be defined using lineage markers; this is then followed by differential expression analysis of additional functional markers. If possible, empirical Bayes methods will again be used to improve power by sharing information on variability across clusters, under the assumption that most functional markers in most clusters are not differentially expressed.



## Statistical approaches

Three main statistical approaches are under development:

- `diffcyt-FDA`: methods based on functional data analysis (FDA), using the [fda](https://cran.r-project.org/web/packages/fda/index.html) package. (Not sure yet whether this can incorporate empirical Bayes.)

- `diffcyt-KS`: methods based on Kolmogorov-Smirnov (KS) tests, and permutation null distributions. We calculate the maximum KS statistic between any two samples, followed by random permutations of group membership labels and re-calculation of the KS statistic to generate a null distribution. (Not sure yet whether this can incorporate empirical Bayes.)

- `diffcyt-med`: simple comparison of medians, i.e. median expression of functional markers within each cluster. (Using empirical Bayes to improve power.)


