#' Test for differential states: method 'diffcyt-DS-LMM'
#' 
#' Calculate tests for differential states within cell populations using method
#' 'diffcyt-DS-LMM'
#' 
#' Calculates tests for differential states within cell populations (i.e. differential
#' expression of cell state markers within clusters), using linear mixed models (LMMs).
#' Clusters are defined using cell type markers, and cell states are characterized by the
#' median transformed expression of cell state markers.
#' 
#' This methodology was originally developed and described by Nowicka et al. (2017),
#' \emph{F1000Research}, and has been modified here to make use of high-resolution
#' clustering to enable investigation of rare cell populations. Note that unlike the
#' original method by Nowicka et al., we do not attempt to manually merge clusters into
#' canonical cell populations. Instead, results are reported at the high-resolution
#' cluster level, and the interpretation of significant differential clusters is left to
#' the user via visualizations such as heatmaps (see the package vignette for an example).
#' 
#' This method fits linear mixed models (LMMs) for each cluster-marker combination (cell
#' state markers only), and calculates differential tests separately for each
#' cluster-marker combination. The response variable in each model is the median
#' arcsinh-transformed marker expression of the cell state marker, which is assumed to
#' follow a Gaussian distribution. There is one model per cluster per cell state marker.
#' Within each model, sample-level weights are included for the number of cells per
#' sample; these weights represent the relative uncertainty in calculating each median
#' value. (Additional uncertainty exists due to variation in the total number of cells per
#' cluster; however, it is not possible to account for this, since there are separate
#' models for each cluster-marker combination.) We also include a filtering step to remove
#' clusters with very small numbers of cells, to improve statistical power.
#' 
#' For more details on the statistical methodology, see Nowicka et al. (2017),
#' \emph{F1000Research} (section 'Differential analysis of marker expression stratified by
#' cell population'.)
#' 
#' The experimental design must be specified using a model formula, which can be created
#' with \code{\link{createFormula}}. Flexible experimental designs are possible, including
#' blocking (e.g. paired designs), batch effects, and continuous covariates. Blocking
#' variables can be included as either random intercept terms or fixed effect terms (see
#' \code{\link{createFormula}}). For paired designs, we recommend using random intercept
#' terms to improve statistical power; see Nowicka et al. (2017), \emph{F1000Research} for
#' details. Batch effects and continuous covariates should be included as fixed effects.
#' 
#' If no random intercept terms are included in the model formula, model fitting is
#' performed using a linear model (LM) instead of a LMM.
#' 
#' The contrast matrix specifying the contrast of interest can be created with
#' \code{\link{createContrast}}. See \code{\link{createContrast}} for more details.
#' 
#' Filtering: Clusters are kept for differential testing if they have at least
#' \code{min_cells} cells in at least \code{min_samples} samples. This removes clusters
#' with very low cell counts across conditions, to improve power.
#' 
#' Weights: Cluster cell counts are used as precision weights within each model (across
#' samples only, i.e. within the model for each cluster); these represent the relative
#' uncertainty in calculating each median value (within each model). See above for
#' details.
#' 
#' 
#' @param d_counts \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'   containing cluster cell counts, from \code{\link{calcCounts}}.
#' 
#' @param d_medians \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'   containing cluster medians (median expression of each marker for each cluster-sample
#'   combination), from \code{\link{calcMedians}}.
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
#'   is appropriate for two-group comparisons (of equal size). Clusters are kept for
#'   differential testing if they have at least \code{min_cells} cells in at least
#'   \code{min_samples} samples.
#' 
#' 
#' @return Returns a new \code{\link[SummarizedExperiment]{SummarizedExperiment}} object,
#'   where rows = cluster-marker combinations, and columns = samples. In the rows,
#'   clusters are repeated for each cell state marker (i.e. the sheets or 'assays' from
#'   the previous \code{d_medians} object are stacked into a single matrix). Differential
#'   test results are stored in the \code{rowData} slot. Results include raw p-values and
#'   adjusted p-values, which can be used to rank cluster-marker combinations by evidence
#'   for differential states within cell populations. The results can be accessed with the
#'   \code{\link[SummarizedExperiment]{rowData}} accessor function.
#' 
#' 
#' @importFrom SummarizedExperiment assay assays rowData 'rowData<-' colData 'colData<-'
#' @importFrom lme4 lmer
#' @importFrom multcomp glht
#' @importFrom stats lm p.adjust
#' @importFrom methods as is
#' 
#' @export
#' 
#' @examples
#' # For a full workflow example demonstrating the use of each function in the 'diffcyt'
#' # pipeline, see the package vignette.
#' 
#' # Create some random data (without differential signal)
#' cofactor <- 5
#' set.seed(123)
#' d_input <- list(
#'   sample1 = sinh(matrix(rnorm(20000, mean = 0, sd = 1), ncol = 20)) * cofactor, 
#'   sample2 = sinh(matrix(rnorm(20000, mean = 0, sd = 1), ncol = 20)) * cofactor, 
#'   sample3 = sinh(matrix(rnorm(20000, mean = 0, sd = 1), ncol = 20)) * cofactor, 
#'   sample4 = sinh(matrix(rnorm(20000, mean = 0, sd = 1), ncol = 20)) * cofactor
#' )
#' # Add differential signal (for some cells and markers in one group)
#' ix_rows <- 901:1000
#' ix_cols <- c(6:10, 16:20)
#' d_input[[3]][ix_rows, ix_cols] <- sinh(matrix(rnorm(1000, mean = 2, sd = 1), ncol = 10)) * cofactor
#' d_input[[4]][ix_rows, ix_cols] <- sinh(matrix(rnorm(1000, mean = 2, sd = 1), ncol = 10)) * cofactor
#' 
#' sample_info <- data.frame(
#'   sample = factor(paste0("sample", 1:4)), 
#'   group = factor(c("group1", "group1", "group2", "group2")), 
#'   stringsAsFactors = FALSE
#' )
#' 
#' marker_info <- data.frame(
#'   marker_name = paste0("marker", 1:20), 
#'   is_marker = rep(TRUE, 20), 
#'   is_type_marker = c(rep(TRUE, 10), rep(FALSE, 10)), 
#'   is_state_marker = c(rep(FALSE, 10), rep(TRUE, 10)), 
#'   stringsAsFactors = FALSE
#' )
#' 
#' # Prepare data
#' d_se <- prepareData(d_input, sample_info, marker_info)
#' # Transform data
#' d_se <- transformData(d_se)
#' # Generate clusters
#' d_se <- generateClusters(d_se)
#' 
#' # Calculate counts
#' d_counts <- calcCounts(d_se)
#' 
#' # Calculate medians (by sample)
#' d_medians <- calcMedians(d_se)
#' 
#' # Create model formula
#' formula <- createFormula(sample_info, cols_fixed = 2, cols_random = 1)
#' # Create contrast matrix
#' contrast <- createContrast(c(0, 1))
#' 
#' # Test for differential states (DS) within clusters
#' res <- testDS_LMM(d_counts, d_medians, formula, contrast)
#' 
testDS_LMM <- function(d_counts, d_medians, formula, contrast, 
                       min_cells = 3, min_samples = NULL) {
  
  if (is.null(min_samples)) {
    min_samples <- ncol(d_counts) / 2
  }
  
  # vector identifying state markers
  id_state_markers <- metadata(d_medians)$id_state_markers
  
  # note: counts are only required for filtering
  counts <- assay(d_counts)
  cluster <- rowData(d_counts)$cluster
  
  # filtering: keep clusters with at least 'min_cells' cells in at least 'min_samples' samples
  tf <- counts >= min_cells
  ix_keep <- apply(tf, 1, function(r) sum(r) >= min_samples)
  
  counts <- counts[ix_keep, , drop = FALSE]
  cluster <- cluster[ix_keep]
  
  # total cell counts per sample (after filtering) (for weights in model fitting)
  n_cells_smp <- colSums(counts)
  
  # LMM/LM testing pipeline
  
  # transpose contrast matrix if created with 'createContrast' (required by 'glht')
  if (ncol(contrast) == 1 & nrow(contrast) > 1) {
    contrast <- t(contrast)
  }
  
  # extract medians and create concatenated matrix
  state_names <- names(assays(d_medians))[id_state_markers]
  meds <- do.call("rbind", {
    lapply(as.list(assays(d_medians)[state_names]), function(a) a[cluster, , drop = FALSE])
  })
  
  meds_all <- do.call("rbind", as.list(assays(d_medians)[state_names]))
  
  # fit models: separate model for each cluster-marker combination
  
  p_vals <- rep(NA, nrow(meds))
  
  for (i in seq_len(nrow(meds))) {
    p_vals[i] <- tryCatch({
      # data for cluster-marker i
      # note: response values are medians
      y <- meds[i, ]
      data_i <- cbind(y, n_cells_smp, formula$data)
      # fit LMM if model formula contains any random effect terms; LM otherwise
      if (formula$random_terms) {
        fit <- lmer(formula$formula, data = data_i, weights = n_cells_smp)
      } else {
        fit <- lm(formula$formula, data = data_i, weights = n_cells_smp)
      }
      # test contrast
      test <- glht(fit, contrast)
      # return p-value
      summary(test)$test$pvalues
      # return NA as p-value if there is an error
    }, error = function(e) NA)
  }
  
  # adjusted p-values (false discovery rates)
  p_adj <- p.adjust(p_vals, method = "fdr")
  
  stopifnot(length(p_vals) == length(p_adj))
  
  # return results in 'rowData' of new 'SummarizedExperiment' object
  
  out <- data.frame(p_vals, p_adj, stringsAsFactors = FALSE)
  
  # fill in any missing rows (filtered clusters) with NAs
  row_data <- as.data.frame(matrix(as.numeric(NA), 
                                   nrow = nlevels(cluster) * length(state_names), 
                                   ncol = ncol(out)))
  colnames(row_data) <- colnames(out)
  
  cluster_nm <- as.numeric(cluster)
  s <- seq(0, nlevels(cluster) * (length(state_names) - 1), by = nlevels(cluster))
  r1 <- rep(cluster_nm, length(state_names))
  r2 <- rep(s, each = length(cluster_nm))
  
  stopifnot(length(s) == length(state_names))
  stopifnot(length(r1) == length(r2))
  
  rows <- r1 + r2
  row_data[rows, ] <- out
  
  # include cluster IDs and marker names
  clus <- factor(rep(levels(cluster), length(state_names)), levels = levels(cluster))
  stat <- factor(rep(state_names, each = length(levels(cluster))), levels = state_names)
  stopifnot(length(clus) == nrow(row_data), length(stat) == nrow(row_data))
  
  row_data <- cbind(data.frame(cluster = clus, marker = stat, stringsAsFactors = FALSE), 
                    row_data)
  
  col_data <- colData(d_medians)
  
  # return object
  res <- SummarizedExperiment(
    meds_all, 
    rowData = row_data, 
    colData = col_data
  )
  
  res
}


