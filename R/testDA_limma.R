#' Test for differential abundance: method 'diffcyt-DA-limma'
#' 
#' Calculate tests for differential abundance of cell populations using method
#' 'diffcyt-DA-limma'
#' 
#' Calculates tests for differential abundance of clusters, using functions from the
#' \code{\link[limma]{limma}} package.
#' 
#' This method uses the \code{\link[limma]{limma}} package (Ritchie et al. 2015,
#' \emph{Nucleic Acids Research}) to fit linear models and calculate empirical Bayes
#' moderated tests at the cluster level. Empirical Bayes methods improve statistical power
#' by sharing information on variability (i.e. variance across samples for a single
#' cluster) between clusters. Since count data are often heteroscedastic, we use the
#' \code{\link[limma]{voom}} method (Law et al. 2014, \emph{Genome Biology}) to transform
#' the raw cluster cell counts and estimate observation-level weights to stabilize the
#' mean-variance relationship. Diagnostic plots are shown if \code{plot = TRUE}.
#' 
#' The experimental design must be specified using a design matrix, which can be created
#' with \code{\link{createDesignMatrix}}. Flexible experimental designs are possible,
#' including blocking (e.g. paired designs), batch effects, and continuous covariates. See
#' \code{\link{createDesignMatrix}} for more details.
#' 
#' For paired designs, either fixed effects or random effects can be used. Fixed effects
#' are simpler, but random effects may improve power in data sets with unbalanced designs
#' or very large numbers of samples. To use fixed effects, provide the block IDs (e.g.
#' patient IDs) to \code{\link{createDesignMatrix}}. To use random effects, provide the
#' \code{block_IDs} argument here instead. This will make use of the \code{limma}
#' \code{\link[limma]{duplicateCorrelation}} methodology. Note that multiple measures per
#' sample are not possible in this case (fixed effects should be used instead). Block IDs
#' should not be included in the design matrix if the \code{limma}
#' \code{\link[limma]{duplicateCorrelation}} methodology is used.
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
#' @param block_IDs (Optional) Vector or factor of block IDs (e.g. patient IDs) for paired
#'   experimental designs, to be included as random effects. If provided, the block IDs
#'   will be included as random effects using the \code{limma}
#'   \code{\link[limma]{duplicateCorrelation}} methodology. Alternatively, block IDs can
#'   be included as fixed effects in the design matrix (\code{\link{createDesignMatrix}}).
#'   See details.
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
#' @param plot Whether to save diagnostic plots for the \code{limma}
#'   \code{\link[limma]{voom}} transformations. Default = TRUE.
#' 
#' @param path Path for diagnostic plots. Default = current working directory.
#' 
#' 
#' @return Returns a new \code{\link[SummarizedExperiment]{SummarizedExperiment}} object,
#'   with differential test results stored in the \code{rowData} slot. Results include raw
#'   p-values and adjusted p-values from the \code{limma} empirical Bayes moderated tests,
#'   which can be used to rank clusters by evidence for differential abundance. The
#'   results can be accessed with the \code{\link[SummarizedExperiment]{rowData}} accessor
#'   function.
#' 
#' 
#' @importFrom SummarizedExperiment assays rowData 'rowData<-' colData 'colData<-'
#' @importFrom limma contrasts.fit voom duplicateCorrelation lmFit eBayes plotSA topTable
#' @importFrom edgeR calcNormFactors DGEList
#' @importFrom methods as is
#' @importFrom grDevices pdf
#' @importFrom graphics plot
#' 
#' @export
#' 
#' @examples
#' # A full workflow example demonstrating the use of each function in the 'diffcyt'
#' # pipeline on an experimental data set is available in the package vignette.
#' 
testDA_limma <- function(d_counts, design, contrast, block_IDs = NULL, 
                         min_cells = 3, min_samples = NULL, 
                         normalize = FALSE, norm_factors = "TMM", 
                         plot = TRUE, path = ".") {
  
  if (!is.null(block_IDs) & !is.factor(block_IDs)) {
    block_IDs <- factor(block_IDs, levels = unique(block_IDs))
  }
  
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
  
  # limma-voom pipeline
  
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
  
  # voom transformation and weights
  if (plot) pdf(file.path(path, "voom_before.pdf"), width = 6, height = 6)
  v <- voom(y, design, plot = TRUE)
  if (plot) dev.off()
  
  # estimate correlation between paired samples
  # (note: paired designs only; >2 measures per sample not allowed)
  if (!is.null(block_IDs)) {
    dupcor <- duplicateCorrelation(v, design, block = block_IDs)
  }
  
  # fit linear models
  if (!is.null(block_IDs)) {
    message("Fitting linear models with random effects term for 'block_IDs'.")
    vfit <- lmFit(v, design, block = block_IDs, correlation = dupcor$consensus.correlation)
  } else {
    vfit <- lmFit(v, design)
  }
  vfit <- contrasts.fit(vfit, contrast)
  
  # calculate empirical Bayes moderated tests
  efit <- eBayes(vfit)
  
  if (plot) pdf(file.path(path, "voom_after.pdf"), width = 6, height = 6)
  plotSA(efit)
  if (plot) dev.off()
  
  # results
  top <- topTable(efit, coef = 1, number = Inf, adjust.method = "BH", sort.by = "none")
  
  if (!all(rownames(top) == cluster)) {
    stop("cluster labels do not match")
  }
  
  # return results in 'rowData' of new 'SummarizedExperiment' object
  
  # fill in any missing rows (filtered clusters) with NAs
  row_data <- as.data.frame(matrix(as.numeric(NA), nrow = nlevels(cluster), ncol = ncol(top)))
  colnames(row_data) <- colnames(top)
  cluster_nm <- as.numeric(cluster)
  row_data[cluster_nm, ] <- top
  
  row_data <- cbind(data.frame(cluster = as.numeric(levels(cluster))), row_data)
  
  res <- d_counts
  
  rowData(res) <- row_data
  
  res
}


