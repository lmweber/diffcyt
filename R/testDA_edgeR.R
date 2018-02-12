#' Test for differential abundance: method 'diffcyt-DA-edgeR'
#' 
#' Calculate tests for differential abundance of cell populations using method
#' 'diffcyt-DA-edgeR'
#' 
#' Calculates tests for differential abundance of clusters, using functions from the
#' \code{\link[edgeR]{edgeR}} package.
#' 
#' This method uses the \code{\link[edgeR]{edgeR}} package (Robinson et al. 2010,
#' \emph{Bioinformatics}; McCarthy et al. 2012, \emph{Nucleic Acids Research}) to fit
#' linear models and calculate empirical Bayes moderated tests at the cluster level.
#' Empirical Bayes methods improve statistical power by sharing information on variability
#' (i.e. variance across samples for a single cluster) between clusters. The statistical
#' methods implemented in the \code{edgeR} package were originally designed for the
#' analysis of gene expression data such as RNA-sequencing counts. Here, we apply these
#' methods to cluster cell counts.
#' 
#' The experimental design must be specified using a design matrix, which can be created
#' with \code{\link{createDesignMatrix}}. Flexible experimental designs are possible,
#' including blocking (e.g. paired designs), batch effects, and continuous covariates. See
#' \code{\link{createDesignMatrix}} for more details.
#' 
#' The contrast matrix specifying the contrast of interest can be created with
#' \code{\link{createContrast}}. See \code{\link{createContrast}} for more details.
#' 
#' Filtering: Clusters are kept for differential testing if they have at least
#' \code{min_cells} cells in at least \code{min_samples} samples. This removes clusters
#' with very low cell counts across conditions, which improves power.
#' 
#' Normalization: Optional normalization factors can be included to adjust for composition
#' effects in the total number of counts per sample (library sizes). For example, if one
#' cell population is more abundant in one condition, while all other populations have
#' equal numbers of cells across conditions, then this effectively reduces the
#' \emph{relative} abundance of the unchanged populations in the first condition, creating
#' false positive differential abundance signals for these populations. Normalization
#' factors can be provided directly as a vector of values representing relative total
#' abundances per sample (where values >1 indicate extra cells in a sample; note this is
#' the inverse of normalization factors as defined in the \code{edgeR} package).
#' Alternatively, normalization factors can be calculated automatically using the 'trimmed
#' mean of M-values' (TMM) method from the \code{edgeR} package (Robinson and Oshlack,
#' 2010). The TMM method assumes that most populations are not differentially abundant;
#' see the \code{edgeR} User's Guide for more details.
#' 
#' 
#' @param d_counts \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'   containing cluster cell counts, from \code{\link{calcCounts}}.
#' 
#' @param design Design matrix, created with \code{\link{createDesignMatrix}}. See
#'   \code{\link{createDesignMatrix}} for details.
#' 
#' @param contrast Contrast matrix, created with \code{\link{createContrast}}. See
#'   \code{\link{createContrast}} for details.
#' 
#' @param min_cells Filtering parameter. Default = 3. Clusters are kept for differential
#'   testing if they have at least \code{min_cells} cells in at least \code{min_samples}
#'   samples.
#' 
#' @param min_samples Filtering parameter. Default = \code{number of samples / 2}, which
#'   is appropriate for two-group comparisons. Clusters are kept for differential testing
#'   if they have at least \code{min_cells} cells in at least \code{min_samples} samples.
#' 
#' @param normalize Whether to include optional normalization factors to adjust for
#'   composition effects (see details). Default = FALSE.
#' 
#' @param norm_factors Normalization factors to use, if \code{normalize = TRUE}. Default =
#'   \code{"TMM"}, in which case normalization factors are calculated automatically using
#'   the 'trimmed mean of M-values' (TMM) method from the \code{edgeR} package.
#' 
#' 
#' @return Returns a new \code{\link[SummarizedExperiment]{SummarizedExperiment}} object,
#'   with differential test results stored in the \code{rowData} slot. Results include raw
#'   p-values and adjusted p-values from the \code{edgeR} empirical Bayes moderated tests,
#'   which can be used to rank clusters by evidence for differential abundance. The
#'   results can be accessed with the \code{\link[SummarizedExperiment]{rowData}} accessor
#'   function.
#' 
#' 
#' @importFrom SummarizedExperiment assays rowData 'rowData<-' colData 'colData<-'
#' @importFrom edgeR calcNormFactors DGEList estimateDisp glmFit glmLRT topTags
#' @importFrom methods as is
#' 
#' @export
#' 
#' @examples
#' # A full workflow example demonstrating the use of each function in the 'diffcyt'
#' # pipeline on an experimental data set is available in the package vignette.
#' 
testDA_edgeR <- function(d_counts, design, contrast, 
                         min_cells = 3, min_samples = NULL, 
                         normalize = FALSE, norm_factors = "TMM") {
  
  if (is.null(min_samples)) {
    min_samples <- ncol(d_counts) / 2
  }
  
  counts <- assays(d_counts)[[1]]
  cluster <- rowData(d_counts)$cluster
  
  # filtering: keep clusters with at least 'min_cells' cells in at least 'min_samples' samples
  tf <- counts >= min_cells
  ix_keep <- apply(tf, 1, function(r) sum(r) >= min_samples)
  
  counts <- counts[ix_keep, ]
  cluster <- cluster[ix_keep]
  
  # edgeR pipeline
  
  # normalization factors
  if (normalize & norm_factors == "TMM") {
    norm_factors <- calcNormFactors(counts, method = "TMM")
  } else if (normalize & norm_factors != "TMM") {
    # edgeR and limma define normalization factors as inverse and require product = 1
    norm_factors <- 1 / norm_factors
    # using geometric mean
    norm_factors <- norm_factors / (prod(norm_factors) ^ (1 / length(norm_factors)))
  }
  
  if (normalize) {
    y <- DGEList(counts, norm.factors = norm_factors)
  } else {
    y <- DGEList(counts)
  }
  
  # estimate dispersions
  # (note: using 'trend.method = "none"')
  y <- estimateDisp(y, design, trend.method = "none")
  
  # fit models
  fit <- glmFit(y, design)
  
  # likelihood ratio tests
  lrt <- glmLRT(fit, contrast = contrast)
  
  # results
  top <- topTags(lrt, n = Inf, adjust.method = "BH", sort.by = "none")
  
  if (!all(rownames(top) == cluster)) {
    stop("cluster labels do not match")
  }
  
  # return results in 'rowData' of new 'SummarizedExperiment' object
  
  # fill in any missing rows (filtered clusters) with NAs
  row_data <- as.data.frame(matrix(as.numeric(NA), nrow = nlevels(cluster), ncol = ncol(top)))
  colnames(row_data) <- colnames(top)
  cluster_nm <- as.numeric(cluster)
  row_data[cluster_nm, ] <- top$table
  
  row_data <- cbind(data.frame(cluster = as.numeric(levels(cluster))), row_data)
  
  res <- d_counts
  
  rowData(res) <- row_data
  
  res
}


