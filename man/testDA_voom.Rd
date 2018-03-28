% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/testDA_voom.R
\name{testDA_voom}
\alias{testDA_voom}
\title{Test for differential abundance: method 'diffcyt-DA-voom'}
\usage{
testDA_voom(d_counts, design, contrast, block_IDs = NULL, min_cells = 3,
  min_samples = NULL, normalize = FALSE, norm_factors = "TMM",
  plot = TRUE, path = ".")
}
\arguments{
\item{d_counts}{\code{\link[SummarizedExperiment]{SummarizedExperiment}} object
containing cluster cell counts, from \code{\link{calcCounts}}.}

\item{design}{Design matrix, created with \code{\link{createDesignMatrix}}. See
\code{\link{createDesignMatrix}} for details.}

\item{contrast}{Contrast matrix, created with \code{\link{createContrast}}. See
\code{\link{createContrast}} for details.}

\item{block_IDs}{(Optional) Vector or factor of block IDs (e.g. patient IDs) for paired
experimental designs, to be included as random effects. If provided, the block IDs
will be included as random effects using the \code{limma}
\code{\link[limma]{duplicateCorrelation}} methodology. Alternatively, block IDs can
be included as fixed effects in the design matrix (\code{\link{createDesignMatrix}}).
See details.}

\item{min_cells}{Filtering parameter. Default = 3. Clusters are kept for differential
testing if they have at least \code{min_cells} cells in at least \code{min_samples}
samples.}

\item{min_samples}{Filtering parameter. Default = \code{number of samples / 2}, which
is appropriate for two-group comparisons (of equal size). Clusters are kept for
differential testing if they have at least \code{min_cells} cells in at least
\code{min_samples} samples.}

\item{normalize}{Whether to include optional normalization factors to adjust for
composition effects (see details). Default = FALSE.}

\item{norm_factors}{Normalization factors to use, if \code{normalize = TRUE}. Default =
\code{"TMM"}, in which case normalization factors are calculated automatically using
the 'trimmed mean of M-values' (TMM) method from the \code{edgeR} package.}

\item{plot}{Whether to save diagnostic plots for the \code{limma}
\code{\link[limma]{voom}} transformations. Default = TRUE.}

\item{path}{Path for diagnostic plots. Default = current working directory.}
}
\value{
Returns a new \code{\link[SummarizedExperiment]{SummarizedExperiment}} object,
  with differential test results stored in the \code{rowData} slot. Results include raw
  p-values and adjusted p-values from the \code{limma} moderated tests, which can be
  used to rank clusters by evidence for differential abundance. The results can be
  accessed with the \code{\link[SummarizedExperiment]{rowData}} accessor function.
}
\description{
Calculate tests for differential abundance of cell populations using method
'diffcyt-DA-voom'
}
\details{
Calculates tests for differential abundance of clusters, using functions from the
\code{\link[limma]{limma}} package and \code{\link[limma]{voom}} method.

This method uses the \code{\link[limma]{limma}} package (Ritchie et al. 2015,
\emph{Nucleic Acids Research}) to fit models and calculate moderated tests at the
cluster level. Moderated tests improve statistical power by sharing information on
variability (i.e. variance across samples for a single cluster) between clusters. Since
count data are often heteroscedastic, we use the \code{\link[limma]{voom}} method (Law
et al. 2014, \emph{Genome Biology}) to transform the raw cluster cell counts and
estimate observation-level weights to stabilize the mean-variance relationship.
Diagnostic plots are shown if \code{plot = TRUE}.

The experimental design must be specified using a design matrix, which can be created
with \code{\link{createDesignMatrix}}. Flexible experimental designs are possible,
including blocking (e.g. paired designs), batch effects, and continuous covariates. See
\code{\link{createDesignMatrix}} for more details.

For paired designs, either fixed effects or random effects can be used. Fixed effects
are simpler, but random effects may improve power in data sets with unbalanced designs
or very large numbers of samples. To use fixed effects, provide the block IDs (e.g.
patient IDs) to \code{\link{createDesignMatrix}}. To use random effects, provide the
\code{block_IDs} argument here instead. This will make use of the \code{limma}
\code{\link[limma]{duplicateCorrelation}} methodology. Note that >2 measures per sample
are not possible in this case (fixed effects should be used instead). Block IDs should
not be included in the design matrix if the \code{limma}
\code{\link[limma]{duplicateCorrelation}} methodology is used.

The contrast matrix specifying the contrast of interest can be created with
\code{\link{createContrast}}. See \code{\link{createContrast}} for more details.

Filtering: Clusters are kept for differential testing if they have at least
\code{min_cells} cells in at least \code{min_samples} samples. This removes clusters
with very low cell counts across conditions, to improve power.

Normalization for the total number of cells per sample (library sizes) and total number
of cells per cluster is automatically performed by the \code{limma} and \code{voom}
functions. Optional normalization factors can also be included to adjust for
composition effects in the total number of counts per sample (library sizes). For
example, if one cell population is more abundant in one condition, while all other
populations have equal numbers of cells across conditions, then this effectively
reduces the \emph{relative} abundance of the unchanged populations in the first
condition, creating false positive differential abundance signals for these
populations. Normalization factors can be provided directly as a vector of values
representing relative total abundances per sample (where values >1 indicate extra cells
in a sample; note this is the inverse of normalization factors as defined in the
\code{edgeR} package). Alternatively, normalization factors can be calculated
automatically using the 'trimmed mean of M-values' (TMM) method from the \code{edgeR}
package (Robinson and Oshlack, 2010). The TMM method assumes that most populations are
not differentially abundant; see the \code{edgeR} User's Guide for more details.
}
\examples{
# For a full workflow example demonstrating the use of each function in the 'diffcyt'
# pipeline, see the package vignette.

# Create some random data (without differential signal)
cofactor <- 5
set.seed(123)
d_input <- list(
  sample1 = sinh(matrix(rnorm(20000, mean = 0, sd = 1), ncol = 20)) * cofactor, 
  sample2 = sinh(matrix(rnorm(20000, mean = 0, sd = 1), ncol = 20)) * cofactor, 
  sample3 = sinh(matrix(rnorm(20000, mean = 0, sd = 1), ncol = 20)) * cofactor, 
  sample4 = sinh(matrix(rnorm(20000, mean = 0, sd = 1), ncol = 20)) * cofactor
)
# Add differential signal (for some cells and markers in one group)
ix_rows <- 901:1000
ix_cols <- 11:20
d_input[[3]][ix_rows, ix_cols] <- sinh(matrix(rnorm(1000, mean = 2, sd = 1), ncol = 10)) * cofactor
d_input[[4]][ix_rows, ix_cols] <- sinh(matrix(rnorm(1000, mean = 2, sd = 1), ncol = 10)) * cofactor

sample_info <- data.frame(
  sample_IDs = paste0("sample", 1:4), 
  group_IDs = factor(c("group1", "group1", "group2", "group2"))
)

marker_info <- data.frame(
  marker_names = paste0("marker", 1:20), 
  is_marker = rep(TRUE, 20), 
  is_type_marker = c(rep(TRUE, 10), rep(FALSE, 10)), 
  is_state_marker = c(rep(FALSE, 10), rep(TRUE, 10))
)

# Prepare data
d_se <- prepareData(d_input, sample_info, marker_info)
# Transform data
d_se <- transformData(d_se)
# Generate clusters
d_se <- generateClusters(d_se)

# Calculate counts
d_counts <- calcCounts(d_se)

# Create design matrix
design <- createDesignMatrix(sample_info, cols_include = 2)
# Create contrast matrix
contrast <- createContrast(c(0, 1))

# Test for differential abundance (DA) of clusters
res <- testDA_voom(d_counts, design, contrast, plot = FALSE)

}