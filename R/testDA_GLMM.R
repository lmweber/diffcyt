#' Test for differential abundance: method 'diffcyt-DA-GLMM'
#' 
#' Calculate tests for differential abundance of clusters using method 'diffcyt-DA-GLMM'
#' 
#' Calculates tests for differential abundance of clusters, using generalized linear mixed
#' models (GLMMs) for each cluster. This methodology was originally developed and
#' described by Nowicka et al. (2017), \emph{F1000Research}.
#' 
#' For more details on the underlying statistical methodology, refer to the paper by
#' Nowicka et al. (2017), \emph{F1000Research}. The implementation here contains several
#' additional modifications. In particular, we use high-resolution clustering and report
#' results at the high-resolution cluster level, instead of relying on a manual
#' cluster-merging step. We also include a filtering step to remove high-resolution
#' clusters with very small numbers of cells, to improve power.
#' 
#' The experimental design must be specified using a model formula, which can be created
#' with \code{\link{createFormula}}. Flexible experimental designs are possible, including
#' batch effects, continuous covariates, and blocking (e.g. for paired designs). Random
#' intercept terms are included for blocks (e.g. paired designs), as well as
#' 'observation-level random effects' to account for overdispersion typically seen in
#' high-dimensional cytometry data (Nowicka et al., 2017, \emph{F1000Research}). See
#' \code{\link{createFormula}} for more details.
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
#' @param formula Model formula object, created with \code{\link{createFormula}}. This
#'   should be a list containing two elements: \code{formula} and \code{data}; the model
#'   formula and data frame of corresponding variables. See \code{\link{createFormula}}
#'   for details.
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
#'   p-values and adjusted p-values, which can be used to rank clusters by evidence for
#'   differential abundance. The results can be accessed with the
#'   \code{\link[SummarizedExperiment]{rowData}} accessor function.
#' 
#' 
#' @importFrom SummarizedExperiment assays rowData 'rowData<-' colData 'colData<-'
#' @importFrom lme4 glmer
#' @importFrom multcomp glht
#' @importFrom IHW ihw adj_pvalues
#' @importFrom methods as is
#' 
#' @export
#' 
#' @seealso to do
#' 
#' @examples
#' # to do
#' 
testDA_GLMM <- function(d_counts, formula, contrast, 
                        min_cells = 3, min_samples = NULL) {
  
  group_IDs <- colData(d_counts)$group
  
  if (is.null(min_samples)) {
    min_samples <- min(table(group_IDs)) - 1
  }
  
  counts <- assays(d_counts)[[1]]
  cluster <- rowData(d_counts)$cluster
  
  # total cell counts per sample (for weights)
  n_cells_smp <- colSums(counts)
  
  # total cell counts per cluster (for hypothesis weighting)
  n_cells <- rowData(d_counts)$n_cells
  
  
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
  n_cells <- n_cells[ix_keep]
  
  
  # transpose contrast matrix if created with 'createContrast' (required by 'glht')
  if (ncol(contrast) == 1 & nrow(contrast) > 1) {
    contrast <- t(contrast)
  }
  
  
  # fit models: separate model for each cluster
  
  p_vals <- rep(NA, length(cluster))
  
  for (i in seq_along(cluster)) {
    # data for cluster i
    # note: divide by total cell counts per sample to enable weights
    y <- counts[i, ] / n_cells_smp
    data_i <- cbind(y, formula$data)
    # fit model
    fit <- glmer(formula$formula, data = data_i, family = "binomial", weights = n_cells_smp)
    # test contrast
    test <- glht(fit, contrast)
    # return p-values
    p_vals[i] <- summary(test)$test$pvalues
  }
  
  
  # calculate adjusted p-values using Independent Hypothesis Weighting (IHW), with total
  # number of cells per cluster as covariate for IHW
  
  ihw_out <- ihw(p_vals, covariates = n_cells, alpha = 0.1)
  p_adj <- adj_pvalues(ihw_out)
  
  
  # return results in 'rowData' of new 'SummarizedExperiment' object
  
  stopifnot(length(p_vals) == length(p_adj))
  out <- data.frame(p_vals, p_adj)
  
  # fill in missing rows (filtered clusters) with NAs
  row_data <- as.data.frame(matrix(NA, nrow = nlevels(cluster), ncol = ncol(out)))
  colnames(row_data) <- colnames(out)
  cluster_num <- as.numeric(cluster)
  row_data[cluster_num, ] <- out
  
  row_data <- cbind(data.frame(cluster = as.numeric(levels(cluster))), row_data)
  
  # store additional sample information in 'colData'
  col_data <- cbind(colData(d_counts), data.frame(group_IDs))
  
  res <- d_counts
  
  rowData(res) <- row_data
  colData(res) <- col_data
  
  res
}



