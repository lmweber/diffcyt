#' Test for differential abundance: method 'diffcyt-DA-edgeR'
#' 
#' Calculate tests for differential abundance of clusters using method 'diffcyt-DA-edgeR'
#' 
#' Calculates tests for differential abundance of clusters, using empirical Bayes
#' moderation of cluster-level variances to improve power.
#' 
#' This method uses the \code{\link[edgeR]{edgeR}} package (Robinson et al. 2010,
#' \emph{Bioinformatics}; McCarthy et al. 2012, \emph{Nucleic Acids Research}) to fit
#' linear models and calculate empirical Bayes moderated tests at the cluster level.
#' Empirical Bayes methods improve statistical power by sharing information on variability
#' (i.e. variance across samples for a single cluster) between clusters. The statistical
#' methods implemented in the \code{edgeR} package were originally designed for the
#' analysis of digital gene expression data such as RNA-sequencing counts. Here, we apply
#' these methods to cluster cell counts.
#' 
#' The experimental design must be specified using a design matrix, which can be created
#' with \code{\link{createDesignMatrix}}. Flexible experimental designs are possible,
#' including batch effects, continuous covariates, and paired designs. See
#' \code{\link{createDesignMatrix}} for more details.
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
#' @param min_cells Filtering parameter. Default = 3. Clusters are kept for differential
#'   testing if they have at least \code{min_cells} cells in at least \code{min_samples}
#'   samples in at least one condition.
#' 
#' @param min_samples Filtering parameter. Default = \code{min(table(group_IDs)) - 1},
#'   i.e. one less than the number of samples in the smallest group. Clusters are kept for
#'   differential testing if they have at least \code{min_cells} cells in at least
#'   \code{min_samples} samples in at least one condition.
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
#' @importFrom edgeR DGEList estimateDisp glmFit glmLRT topTags
#' @importFrom methods as is
#' 
#' @export
#' 
#' @seealso to do
#' 
#' @examples
#' # to do
#' 
testDA_edgeR <- function(d_counts, design, contrast, 
                         min_cells = 3, min_samples = NULL) {
  
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
  
  # edgeR pipeline
  
  # prepare object
  y <- DGEList(counts)
  # normalization factors not required in this context (cluster cell counts)
  #y <- calcNormFactors(y)
  # estimate dispersions
  y <- estimateDisp(y, design)
  ## also try with 'trend.method = "none"'
  
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
  
  # store additional sample information in 'colData'
  col_data <- cbind(colData(d_counts), data.frame(group_IDs))
  
  res <- d_counts
  
  rowData(res) <- row_data
  colData(res) <- col_data
  
  res
}


