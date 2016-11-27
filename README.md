# diffcyt

R package for differential analysis of cell populations in high-dimensional flow cytometry and mass cytometry (CyTOF) data.

Work in progress.


## Contents

Current version implements the following:


### Differential abundance of cell populations

Differential abundance of cell populations detected by automatic clustering using FlowSOM. Abundance is measured as cell population frequencies.


### Differential expression of functional markers

Differential expression of functional markers in cell populations detected by automatic clustering on lineage markers only; i.e. cell populations are detected by clustering on lineage markers only, followed by analysis of differential expression of additional functional markers.

Multiple methods are available for calculating differential expression of functional markers:

- Comparison of medians; i.e. median expression of the functional marker in each sample. Empirical Bayes methods from the `limma` package are used to share information across clusters to improve power, under the assumption that most functional markers in most clusters are not differentially expressed.

- Maximum Kolmogorov-Smirnov tests with permutation null distributions. We calculate the maximum Kolmogorov-Smirnov test statistic between any sample in group 1 compared to any sample in group 2, followed by random permutations of the group membership labels and re-calculation of the test statistic to determine a null distribution, giving an overall p-value. Currently only 2 groups are possible for these tests.

- Functional data analysis approach. Work in progress.


