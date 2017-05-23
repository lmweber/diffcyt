# diffcyt

R package: Statistical methods for differential discovery in high-dimensional flow cytometry and mass cytometry (CyTOF) data.

Under development.



## Types of differential discovery

We consider two main types of differential discovery / analysis:

- *Differential abundance* of cell populations. Cell populations are defined using unsupervised clustering.

- *Differential expression of functional markers* within cell populations. Cell populations are defined using unsupervised clustering on a subset of protein markers (e.g. surface markers in immunology). Additional functional markers are then analyzed for differential expression within the clusters. The two sets of markers (clustering and functional) must be specified by the user.



## Statistical approaches

Several different statistical approaches are under development:

- `diffcyt-DA`: differential abundance; using empirical Bayes methods (from the `limma` package) to share information on variability between clusters

- `diffcyt-med`: differential expression of functional markers; characterizing each functional marker expression profile by its median, and using empirical Bayes methods to share information on variability between clusters

- `diffcyt-FDA`: differential expression of functional markers; using methods from functional data analysis (FDA) to characterize each functional marker expression profile by its empirical cumulative distribution function (ECDF), which contains significantly more information than a single value (e.g. median). Also using number of cells per cluster-sample combination as weights to represent the uncertainty in calculating each ECDF. Using permutation tests to test for differences between the groups of ECDF curves.

- `diffcyt-KS`: under development; possibly using Kolmogorov-Smirnov (KS) tests to compare the groups of ECDF curves. Permutation null distributions. Calculate the maximum KS statistic between any two samples; random permutations of group membership labels to calculate null distribution.

- `diffcyt-LM`: under development


