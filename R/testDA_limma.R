#' Test for differential abundance: method 'diffcyt-DA-limma'
#' 
#' Calculate tests for differential abundance of clusters using method 'diffcyt-DA-limma'
#' 
#' Calculates tests for differential abundance of clusters, using empirical Bayes
#' moderation of cluster-level variances to improve power.
#' 
#' The \code{\link[limma]{limma}} package (Ritchie et al. 2015, \emph{Nucleic Acids
#' Research}) is used to fit linear models and calculate empirical Bayes moderated tests.
#' Empirical Bayes methods improve statistical power by sharing information on variability
#' (i.e. variance across samples for a single cluster) between clusters. Since count data
#' are often heteroscedastic, we use the  \code{\link[limma]{voom}} method (Law et al.
#' 2014, \emph{Genome Biology}) to transform the raw cluster cell counts and estimate
#' observation-level weights to stabilize the mean-variance relationship. Diagnostic plots
#' are shown if \code{plot = TRUE}.
#' 
#' The experimental design must be specified using a design matrix, which can be created
#' with \code{\link{createDesignMatrix}}. Flexible experimental designs are possible,
#' including batch effects, continuous covariates, and paired designs. See
#' \code{\link{createDesignMatrix}} for more details.
#' 
#' For paired designs, either fixed effects or random effects can be used. Fixed effects
#' are simpler, but random effects may improve power for unbalanced designs. To use fixed
#' effects, provide the block IDs (e.g. patient IDs) to \code{\link{createDesignMatrix}}.
#' To use random effects, provide the \code{block_IDs} argument here instead. This will
#' make use of the \code{limma} \code{\link[limma]{duplicateCorrelation}} methodology;
#' note that block IDs should not be included in the design matrix in this case.
#' 
#' The contrast matrix specifying the contrast of interest can be created with
#' \code{\link{createContrast}}. See \code{\link{createContrast}} for more details.
#' 
#' Filtering: Clusters are kept for differential testing if they have at least
#' \code{min_cells} cells in at least \code{min_samples} samples in at least one
#' condition. This removes clusters with very low cell counts across conditions, which
#' improves power.
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
#'   samples in at least one condition.
#' 
#' @param min_samples Filtering parameter. Default = \code{min(table(group_IDs)) - 1},
#'   i.e. one less than the number of samples in the smallest group. Clusters are kept for
#'   differential testing if they have at least \code{min_cells} cells in at least
#'   \code{min_samples} samples in at least one condition.
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
testDA_limma <- function(d_counts, design, contrast, 
                         block_IDs = NULL, 
                         min_cells = 3, min_samples = NULL, 
                         plot = TRUE, path = ".") {
  
  if (!is.null(block_IDs) & !is.factor(block_IDs)) {
    block_IDs <- factor(block_IDs, levels = unique(block_IDs))
  }
  
  group_IDs <- colData(d_counts)$group
  
  if (is.null(min_samples)) {
    min_samples <- min(table(group_IDs)) - 1
  }
  
  counts <- assays(d_counts)[[1]]
  cluster <- rowData(d_counts)$cluster
  
  # filtering: keep clusters with at least 'min_cells' cells in at least 'min_samples'
  # samples in at least one condition
  ix_keep <- rep(FALSE, length(cluster))
  tf <- counts >= min_cells
  for (g in seq_along(levels(group_IDs))) {
    grp <- group_IDs == levels(group_IDs)[g]
    ix_keep[rowSums(tf[, grp, drop = FALSE]) >= min_samples] <- TRUE
  }
  
  counts <- counts[ix_keep, ]
  cluster <- cluster[ix_keep]
  
  # voom transformation and weights
  if (plot) pdf(file.path(path, "voom_before.pdf"), width = 6, height = 6)
  v <- voom(counts, design, plot = TRUE)
  if (plot) dev.off()
  
  # estimate correlation between paired samples
  # (note: paired designs only; >2 measurements per sample not allowed)
  if (!is.null(block_IDs)) {
    dupcor <- duplicateCorrelation(v, design, block = block_IDs)
  }
  
  # fit linear models
  if (!is.null(block_IDs)) {
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
  
  # store additional sample information in 'colData'
  if (!is.null(block_IDs)) {
    col_data <- cbind(colData(d_counts), block = data.frame(block_IDs))
  }
  
  res <- d_counts
  
  rowData(res) <- row_data
  colData(res) <- col_data
  
  res
}


