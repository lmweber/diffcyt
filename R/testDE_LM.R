#' Test for differential expression within clusters (method: 'diffcyt-LM')
#' 
#' Calculate tests for differential expression of functional markers within clusters 
#' (method: 'diffcyt-LM')
#' 
#' The 'diffcyt-LM' methodology fits a single linear model for each cluster, to test for
#' differences between the ECDF curves for each group of samples (e.g. diseased vs.
#' healthy).
#' 
#' In the case of paired tests, a linear mixed model is used, with a random intercept term
#' for each block (e.g. patient IDs).
#' 
#' By default, the ECDF curves are 'clipped' prior to testing, so that only the more 
#' informative indices in the middle of the curves are used for testing (by default, the
#' middle two-thirds are used).
#' 
#' By using the ECDFs over a range of values, these tests make more use of information
#' than tests that are based on only a single summary value for each cluster-marker
#' combination (e.g. the median).
#' 
#' By default, the numbers of cells per cluster-sample combination are used as weights in
#' the linear (or linear mixed) models. The final p-values by cluster are also adjusted
#' using the Independent Hypothesis Weighting (IHW) method (Ignatiadis et al. 2016) using
#' the total number of cells per cluster as covariate.
#' 
#' Filtering: Clusters are kept for testing if there are at least \code{min_cells} cells 
#' per sample in at least \code{min_samples} samples in either condition. Filtered 
#' clusters are removed from differential expression testing for all markers. Clusters
#' containing any zero values (zero cells per cluster-sample combination) are also
#' removed.
#' 
#' Currently, only two-group comparisons are possible. Future work will extend this method
#' to allow more complex comparisons (contrasts).
#' 
#' Alternative methodologies for testing for differential expression within clusters are 
#' available in the functions \code{\link{testDE_med}}, \code{\link{testDE_KS}}, and 
#' \code{\link{testDE_FDA}}.
#' 
#' 
#' @param d_counts \code{\link[SummarizedExperiment]{SummarizedExperiment}} object 
#'   containing cluster cell counts, from \code{\link{calcCounts}}.
#' 
#' @param d_medians \code{\link[SummarizedExperiment]{SummarizedExperiment}} object 
#'   containing cluster medians (median expression of functional markers), from 
#'   \code{\link{calcMedians}}.
#' 
#' @param d_ecdfs \code{\link[SummarizedExperiment]{SummarizedExperiment}} object 
#'   containing empirical cumulative distribution functions (ECDFs), calculated with 
#'   \code{\link{calcECDFs}}.
#' 
#' @param group_IDs Vector or factor of group membership IDs for each sample (e.g. 
#'   diseased vs. healthy, or treated vs. untreated). Vectors are converted to factors 
#'   internally. The first level of the factor will be used as the reference level for 
#'   differential testing. Currently, only two-group comparisons are implemented.
#' 
#' @param weighted Whether to include weights for the number of cells per cluster-sample
#'   combination. Weights (number of cells) represent the relative uncertainty in 
#'   calculating each ECDF curve. If \code{weighted = FALSE}, unweighted tests will be 
#'   calculated instead, which ignores this source of uncertainty. Note these are relative
#'   weights within each cluster only (not across clusters). Default = TRUE.
#' 
#' @param paired Whether to perform paired tests. Set to TRUE and provide the 
#'   \code{block_IDs} argument (e.g. patient IDs) to calculate paired tests. Default = 
#'   FALSE.
#' 
#' @param block_IDs Vector or factor of block IDs for samples (e.g. patient ID), required 
#'   for paired tests. Default = NULL.
#' 
#' @param min_cells Filtering parameter. Default = 5. Clusters are kept if there are at 
#'   least \code{min_cells} cells per sample in at least \code{min_samples} samples in
#'   either condition. Filtered clusters are removed from differential expression testing
#'   for all markers.
#' 
#' @param min_samples Filtering parameter. Default = \code{n - 1}, where \code{n} = number
#'   of replicates in smallest group. Clusters are kept if there are at least
#'   \code{min_cells} cells per sample in at least \code{min_samples} samples in either
#'   condition. Filtered clusters are removed from differential expression testing for all
#'   markers.
#' 
#' @param clip_ecdfs_low Proportion of ECDF indices (i.e. values at which the ECDF curves
#'   are evaluated) to remove from lower end. The indices remaining between the two
#'   clipped ends (i.e. the more informative values in the middle) are used for testing. 
#'   Default: \code{clip_ecdfs_low = 1/6} and \code{clip_ecdfs_high = 1/6}; middle
#'   two-thirds are kept.
#' 
#' @param clip_ecdfs_high Proportion of ECDF indices (i.e. values at which the ECDF curves
#'   are evaluated) to remove from upper end. The indices remaining between the two 
#'   clipped ends (i.e. the more informative values in the middle) are used for testing. 
#'   Default: \code{clip_ecdfs_low = 1/6} and \code{clip_ecdfs_high = 1/6}; middle 
#'   two-thirds are kept.
#' 
#' 
#' @return Returns new \code{\link[SummarizedExperiment]{SummarizedExperiment}} object, 
#'   where rows = cluster-marker combinations, columns = samples, values = median marker 
#'   expression. In the rows, clusters are repeated for each functional marker (i.e. the 
#'   sheets or 'assays' from the previous \code{d_medians} object are stacked into a 
#'   single matrix). Differential test results are stored in the 'rowData' slot. Results 
#'   include p-values and adjusted p-values (using Independent Hypothesis Weighting; 
#'   Ignatiadis et al. 2016), which can be used to rank cluster-marker combinations by 
#'   evidence for differential expression. The results can be accessed with the 
#'   \code{\link[SummarizedExperiment]{rowData}} accessor function.
#' 
#' 
#' @importFrom SummarizedExperiment assays rowData colData
#' @importFrom stats model.matrix sd lm
#' @importFrom lme4 lmer
#' @importFrom multcomp glht
#' @importFrom IHW ihw adj_pvalues
#' @importFrom reshape2 melt
#' 
#' @export
#' 
#' @examples
#' library(flowCore)
#' library(SummarizedExperiment)
#' 
#' # filenames
#' files <- list.files(system.file("extdata", package = "diffcyt"), 
#'                     pattern = "\\.fcs$", full.names = TRUE)
#' files_BCRXL <- files[grep("BCRXL", files)]
#' files_ref <- files[grep("ref", files)]
#' 
#' # load data
#' files_load <- c(files_BCRXL, files_ref)
#' d_input <- lapply(files_load, read.FCS, transformation = FALSE, truncate_max_range = FALSE)
#' 
#' # sample IDs and group IDs
#' sample_IDs <- gsub("\\.fcs$", "", basename(files_load))
#' sample_IDs
#' group_IDs <- gsub("^patient[0-9]_", "", sample_IDs)
#' group_IDs
#' 
#' # set group reference level for differential testing
#' group_IDs <- factor(group_IDs, levels = c("ref", "BCRXL"))
#' group_IDs
#' 
#' # indices of all marker columns, lineage markers, and functional markers
#' # (see Table 1 in Bruggner et al. 2014)
#' cols_markers <- c(3:4, 7:9, 11:19, 21:22, 24:26, 28:31, 33)
#' cols_lineage <- c(3:4, 9, 11, 12, 14, 21, 29, 31, 33)
#' cols_func <- setdiff(cols_markers, cols_lineage)
#' 
#' # prepare data
#' # (note: using lineage markers for clustering, and functional markers for DE testing)
#' d_se <- prepareData(d_input, sample_IDs, cols_markers, cols_lineage, cols_func)
#' 
#' # transform data
#' d_se <- transformData(d_se, cofactor = 5)
#' 
#' # generate clusters
#' # (note: using small number of clusters for demonstration purposes in this example)
#' d_se <- generateClusters(d_se, xdim = 4, ydim = 4, seed = 123)
#' 
#' # calculate cluster cell counts
#' d_counts <- calcCounts(d_se)
#' 
#' # calculate cluster medians
#' d_medians <- calcMedians(d_se)
#' 
#' # calculate ECDFs
#' d_ecdfs <- calcECDFs(d_se)
#' 
#' # subset marker expresison values
#' d_vals <- subsetVals(d_se)
#' 
#' 
#' #############################################################################
#' # Test for differential expression (DE) of functional markers within clusters
#' # (method 'diffcyt-LM')
#' #############################################################################
#' 
#' # create block IDs for paired tests (this is a paired data set, so we use 1 block per patient)
#' patient_IDs <- factor(gsub("_(BCRXL|ref)$", "", sample_IDs))
#' patient_IDs <- as.numeric(patient_IDs)
#' patient_IDs
#' 
#' # test for differential expression (DE) of functional markers within clusters
#' res_DE <- testDE_LM(d_counts, d_medians, d_ecdfs, group_IDs, 
#'                     paired = TRUE, block_IDs = patient_IDs)
#' 
#' # show results using 'rowData' accessor function
#' rowData(res_DE)
#' 
#' # sort to show top (most highly significant) cluster-marker combinations first
#' head(rowData(res_DE)[order(rowData(res_DE)$p_adj), ], 10)
#' 
testDE_LM <- function(d_counts, d_medians, d_ecdfs, group_IDs, 
                      weighted = TRUE, paired = FALSE, block_IDs = NULL, 
                      min_cells = 5, min_samples = NULL, 
                      clip_ecdfs_low = 1/6, clip_ecdfs_high = 1/6) {
  
  if (!is.factor(group_IDs)) group_IDs <- factor(group_IDs, levels = unique(group_IDs))
  
  if (paired & is.null(block_IDs)) {
    stop("'block_IDs' argument is required for paired tests")
  }
  
  if (is.null(min_samples)) {
    min_samples <- min(table(group_IDs)) - 1
  }
  
  counts <- assays(d_counts)[[1]]
  cluster <- rowData(d_counts)$cluster
  
  func_names <- names(assays(d_ecdfs))
  
  # filtering
  grp_ref <- group_IDs == levels(group_IDs)[1]
  tf <- counts >= min_cells
  ix_keep <- (rowSums(tf[, grp_ref]) >= min_samples) & (rowSums(tf[, !grp_ref]) >= min_samples)
  
  counts <- counts[ix_keep, ]
  cluster <- cluster[ix_keep]
  
  # remove any remaining rows (clusters) with zeros
  ix_zeros <- apply(counts, 1, function(r) any(r == 0))
  
  counts <- counts[!ix_zeros, ]
  cluster <- cluster[!ix_zeros]
  
  # note: assumes same resolution for all ECDFs
  resolution <- length(assays(d_ecdfs)[[1]][1, 1][[1]])
  index <- seq_len(resolution)
  # clip ECDFs (remove lower and upper values)
  index <- index[-c(1:(resolution * clip_ecdfs_low), resolution:((1 - clip_ecdfs_high) * resolution))]
  
  # create design matrix
  index_mm <- rep(factor(index), length(group_IDs))
  grp_mm <- rep(group_IDs, each = length(index))
  X <- model.matrix(~ index_mm + grp_mm)
  
  # set up matrix for p-values
  p_vals <- matrix(NA, nrow = length(levels(cluster)), ncol = length(func_names))
  rownames(p_vals) <- levels(cluster)
  colnames(p_vals) <- func_names
  
  # fit models and calculate p-values (for cluster 'i' and marker 'j')
  for (j in seq_along(func_names)) {
    for (i in cluster) {  # note: only selects from clusters 'i' that passed filter
      
      Y <- assays(d_ecdfs)[[j]][i, ]
      
      # clip ECDFs (remove lower and upper values)
      Y <- lapply(Y, function(y) y[index])
      
      # transformation: mean = 0, sd = 1 (for each index)
      Y <- do.call(cbind, Y)
      rownames(Y) <- index
      Y <- t(apply(Y, 1, function(y) (y - mean(y)) / sd(y)))
      
      Y <- unlist(as.data.frame(Y))
      
      # weights: number of cells per sample
      if (weighted) {
        weights <- rep(counts[i, ], each = length(index))
      } else {
        weights <- NULL
      }
      
      # fit model
      if (paired) {
        block_mm <- factor(rep(block_IDs, each = length(index)))
        fit <- lmer(Y ~ X + (1 | block_mm), weights = weights)
        
        # set up contrast: test significance of group membership term
        contr <- c(rep(0, ncol(X) - 1), 1)
        contr <- matrix(contr, 1)
        
        test <- glht(fit, linfct = contr)
        
        # extract p-value for group membership term
        p_val <- summary(test)$test$pvalues
        
      } else {
        fit <- lm(Y ~ X, weights = weights)
        
        # extract p-value for group membership term (note: assumes this is the final term)
        p_val <- summary(fit)$coefficients[, 4][length(summary(fit)$coefficients[, 4])]
      }
      
      p_vals[i, j] <- p_val
    }
  }
  
  # return results
  res <- as.data.frame(p_vals)
  res$cluster <- levels(cluster)
  
  # calculate adjusted p-values using Independent Hypothesis Weighting (IHW)
  # using total number of cells per cluster as the covariate for IHW
  res$n_cells <- rowData(d_counts)$n_cells
  
  res <- melt(res, id.vars = c("cluster", "n_cells"), 
              variable.name = "marker", 
              value.name = "p_val")
  
  ihw_out <- ihw(p_val ~ n_cells, data = res, alpha = 0.1)
  
  res$p_adj <- adj_pvalues(ihw_out)
  
  # return new 'SummarizedExperiment' object with results stored in 'rowData'
  meds <- do.call("rbind", {
    lapply(as.list(assays(d_medians)[func_names]), function(a) a[levels(cluster), ])
  })
  
  if (!all(rownames(meds) == res$cluster)) {
    stop("cluster labels do not match")
  }
  
  if (paired) {
    col_data <- cbind(colData(d_medians), data.frame(group_IDs), data.frame(block_IDs))
  } else {
    col_data <- cbind(colData(d_medians), data.frame(group_IDs))
  }
  
  row_data <- res
  
  res_DE <- SummarizedExperiment(meds, rowData = row_data, colData = col_data)
  
  res_DE
}


