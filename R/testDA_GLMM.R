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
#' Filtering: Clusters are kept for differential testing if they have at least
#' \code{min_cells} cells in at least \code{min_samples} samples in at least one
#' condition. This removes clusters with very low cell counts across conditions, which
#' improves statistical power.
#' 
#' 
#' @param d_counts \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'   containing cluster cell counts, from \code{\link{calcCounts}}.
#' 
#' @param group_IDs Vector or factor of group membership labels for each sample (e.g.
#'   diseased vs. healthy, or treated vs. untreated). Vectors are converted to factors
#'   internally. The user must ensure that this is identical to the \code{group_IDs}
#'   provided previously to \code{\link{prepareData}}.
#' 
#' @param contrast Contrast specifying the differential comparison of interest for
#'   testing. Each contrast should be in a separate row. If not provided, the default is
#'   to compare the second vs. first level of \code{group_IDs}.
#' 
#' @param batch_IDs (Optional) Vector or factor of batch IDs. Batch effects are removed by
#'   adding the batch IDs to the GLMM.
#' 
#' @param covariates (Optional) Numeric matrix of continuous covariates, to be added to
#'   the GLMM in order to remove their effects.
#' 
#' @param block_IDs (Optional) Vector or factor of block IDs, for paired experimental
#'   designs.
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
#' @importFrom IHW ihw
#' @importFrom stats formula
#' @importFrom methods as is
#' 
#' @export
#' 
#' @seealso to do
#' 
#' @examples
#' # to do
#' 
testDA_GLMM <- function(d_counts, group_IDs, contrast = NULL, 
                        batch_IDs = NULL, covariates = NULL, block_IDs = NULL, 
                        min_cells = 3, min_samples = NULL) {
  
  if (!is.factor(group_IDs)) {
    group_IDs <- factor(group_IDs, levels = unique(group_IDs))
  }
  if (!is.null(batch_IDs) & !is.factor(batch_IDs)) {
    batch_IDs <- factor(batch_IDs, levels = unique(batch_IDs))
  }
  if (!is.null(block_IDs) & !is.factor(block_IDs)) {
    block_IDs <- factor(block_IDs, levels = unique(block_IDs))
  }
  
  if (!is.null(covariates) & !is.matrix(covariates) & !is.numeric(covariates)) {
    stop("'covariates' must be provided as a numeric matrix, with one column for each covariate")
  }
  
  if (is.null(min_samples)) {
    min_samples <- min(table(group_IDs)) - 1
  }
  
  counts <- assays(d_counts)[[1]]
  cluster <- rowData(d_counts)$cluster
  
  sample_IDs <- colData(d_counts)$sample
  stopifnot(all(sample_IDs == colnames(counts)))
  
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
  
  
  # set up model formula
  # (note: allows batch effects and covariates)
  if (!is.null(batch_IDs) & !is.null(covariates)) {
    formula <- y / n_cells_smp ~ group_IDs + (1 | sample_IDs) + (1 | block_IDs) + batch_IDs + covariates
    ### to do: check if covariates work correctly when provided as numeric matrix
  } else if (!is.null(batch_IDs)) {
    formula <- y / n_cells_smp ~ group_IDs + (1 | sample_IDs) + (1 | block_IDs) + batch_IDs
  } else if (!is.null(covariates)) {
    formula <- y / n_cells_smp ~ group_IDs + (1 | sample_IDs) + (1 | block_IDs) + covariates
  } else {
    formula <- y / n_cells_smp ~ group_IDs + (1 | sample_IDs) + (1 | block_IDs)
  }
  
  
  # set up contrasts
  # (note: if not specified, default is to compare 2nd vs. 1st level of 'group_IDs')
  if (is.null(contrast)) {
    levs <- paste0("group_IDs", levels(group_IDs))
    contr_name <- paste(as.character(levs[2]), "-", as.character(levs[1]))
    contr <- rep(0, length(levs))
    contr[2] <- 1
    contrast <- matrix(contr, nrow = 1)
    rownames(contrast) <- contr_name
    colnames(contrast) <- levs
  }
  
  
  # fit models: separate model for each cluster
  p_vals <- rep(NA, length(cluster))
  
  for (i in seq_along(cluster)) {
    # data for cluster i
    y <- counts[i, ]
    # fit model
    fit <- glmer(formula, family = "binomial", weights = n_cells_smp)
    # test contrast
    test <- glht(fit, contrast)
    # return p-values
    p_vals[i] <- summary(test)$test$pvalues
    cat(".")
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



