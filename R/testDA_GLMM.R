#' Test for differential abundance: method 'diffcyt-DA-GLMM'
#' 
#' Calculate tests for differential abundance of cell populations using method
#' 'diffcyt-DA-GLMM'
#' 
#' Calculates tests for differential abundance of clusters, using generalized linear mixed
#' models (GLMMs).
#' 
#' This methodology was originally developed and described by Nowicka et al. (2017),
#' \emph{F1000Research}, and has been modified here to make use of high-resolution
#' clustering to enable investigation of rare cell populations. Note that unlike the
#' original method by Nowicka et al., we do not attempt to manually merge clusters into
#' canonical cell populations. Instead, results are reported at the high-resolution
#' cluster level, and the interpretation of significant differential clusters is left to
#' the user via visualizations such as heatmaps and tSNE plots (see the package vignette
#' for a detailed example).
#' 
#' This method fits generalized linear mixed models (GLMMs) for each cluster, and
#' calculates differential tests separately for each cluster. The response variables in
#' the models are the cluster cell counts, which are assumed to follow a binomial
#' distribution. There is one model per cluster. We also include a filtering step to
#' remove clusters with very small numbers of cells, to improve statistical power.
#' 
#' For more details on the statistical methodology, see Nowicka et al. (2017),
#' \emph{F1000Research} (section 'Differential cell population abundance'.)
#' 
#' The experimental design must be specified using a model formula, which can be created
#' with \code{\link{createFormula}}. Flexible experimental designs are possible, including
#' blocking (e.g. paired designs), batch effects, and continuous covariates. Blocking
#' variables can be included as either random intercept terms or fixed effect terms (see
#' \code{\link{createFormula}}). For paired designs, we recommend using random intercept
#' terms to improve statistical power; see Nowicka et al. (2017), \emph{F1000Research} for
#' details. Batch effects and continuous covariates should be included as fixed effects.
#' In addition, we include random intercept terms for each sample to account for
#' overdispersion typically seen in high-dimensional cytometry count data. The
#' sample-level random intercept terms are known as 'observation-level random effects'
#' (OLREs); see Nowicka et al. (2017), \emph{F1000Research} for more details.
#' 
#' The contrast matrix specifying the contrast of interest can be created with
#' \code{\link{createContrast}}. See \code{\link{createContrast}} for more details.
#' 
#' Filtering: Clusters are kept for differential testing if they have at least
#' \code{min_cells} cells in at least \code{min_samples} samples. This removes clusters
#' with very low cell counts across conditions, which improves power.
#' 
#' 
#' @param d_counts \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'   containing cluster cell counts, from \code{\link{calcCounts}}.
#' 
#' @param formula Model formula object, created with \code{\link{createFormula}}. This
#'   should be a list containing three elements: \code{formula}, \code{data}, and
#'   \code{random_terms}: the model formula, data frame of corresponding variables, and
#'   variable indicating whether the model formula contains any random effect terms. See
#'   \code{\link{createFormula}} for details.
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
#' @importFrom methods as is
#' @importFrom stats p.adjust
#' 
#' @export
#' 
#' @examples
#' # A full workflow example demonstrating the use of each function in the 'diffcyt'
#' # pipeline on an experimental data set is available in the package vignette.
#' 
testDA_GLMM <- function(d_counts, formula, contrast, 
                        min_cells = 3, min_samples = NULL) {
  
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
  
  # total cell counts per sample (after filtering) (for weights in model fitting)
  n_cells_smp <- colSums(counts)
  
  # GLMM testing pipeline
  
  # transpose contrast matrix if created with 'createContrast' (required by 'glht')
  if (ncol(contrast) == 1 & nrow(contrast) > 1) {
    contrast <- t(contrast)
  }
  
  # fit models: separate model for each cluster
  
  p_vals <- rep(NA, length(cluster))
  
  for (i in seq_along(cluster)) {
    # data for cluster i
    # note: divide by total number of cells per sample (after filtering) to get
    # proportions instead of counts
    y <- counts[i, ] / n_cells_smp
    data_i <- cbind(y, n_cells_smp, formula$data)
    # fit model
    # note: provide proportions (y) together with weights for total number of cells per
    # sample (n_cells_smp); this is equivalent to providing counts
    fit <- glmer(formula$formula, data = data_i, family = "binomial", weights = n_cells_smp)
    # test contrast
    test <- glht(fit, contrast)
    # return p-value
    p_vals[i] <- summary(test)$test$pvalues
  }
  
  # adjusted p-values (false discovery rates)
  p_adj <- p.adjust(p_vals, method = "fdr")
  
  stopifnot(length(p_vals) == length(p_adj))
  
  # return results in 'rowData' of new 'SummarizedExperiment' object
  
  out <- data.frame(p_vals, p_adj)
  
  # fill in any missing rows (filtered clusters) with NAs
  row_data <- as.data.frame(matrix(as.numeric(NA), nrow = nlevels(cluster), ncol = ncol(out)))
  colnames(row_data) <- colnames(out)
  cluster_nm <- as.numeric(cluster)
  row_data[cluster_nm, ] <- out
  
  row_data <- cbind(data.frame(cluster = as.numeric(levels(cluster))), row_data)
  
  res <- d_counts
  
  rowData(res) <- row_data
  
  res
}



