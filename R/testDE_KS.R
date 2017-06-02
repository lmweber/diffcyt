#' Test for differential expression within clusters (method: 'diffcyt-KS')
#' 
#' Calculate tests for differential expression of functional markers within clusters 
#' (method: 'diffcyt-KS')
#' 
#' The 'diffcyt-KS' methodology uses Kolmogorov-Smirnov (KS) tests to compare the 
#' functional marker expression profiles in each cluster between the two groups of samples
#' (e.g. diseased vs. healthy).
#' 
#' KS tests are sensitive to differences in both location and shape of the respective 
#' distributions (marker expression profiles). Therefore, this method makes use of more of
#' the information contained in the expression profiles, compared to testing based on 
#' differences in medians only.
#' 
#' For unpaired tests (default), the maximum KS statistic between any two samples between 
#' the two groups is used as the test statistic. Permutations of group membership labels 
#' are used to generate a null distribution. The minimum p-value in this case is \code{1 /
#' n_perm}, where \code{n_perm} is the number of permutations.
#' 
#' For paired tests, p-values are calculated directly from the two-sample KS tests for 
#' each pair. The p-values for the pairs are combined using Fisher's method for combining 
#' p-values ("Fisher's combined probability test"), giving an overall p-value for each 
#' cluster.
#' 
#' Two levels of filtering are included by default. First, for each cluster-marker 
#' combination, samples with less than \code{min_cells} cells are removed. Second, any 
#' clusters with less than \code{min_samples} samples remaining in either condition are 
#' removed (for this marker), with p-value = NA returned in these cases.
#' 
#' Sample sizes or weights (number of cells per cluster-sample combination) are used 
#' during the calculation of the p-values.
#' 
#' Currently, only two-group comparisons are possible with this method.
#' 
#' The \code{\link[BiocParallel]{BiocParallel}} package is used for parallelized 
#' evaluation on multi-core systems to speed up runtime (for unpaired tests only).
#' 
#' Alternative methodologies for testing for differential expression within clusters are 
#' available in the functions \code{\link{testDE_med}}, \code{\link{testDE_KS}}, and 
#' \code{\link{testDE_LM}}.
#' 
#' 
#' @param d_counts \code{\link[SummarizedExperiment]{SummarizedExperiment}} object 
#'   containing cluster cell counts, from \code{\link{calcCounts}}. Required for filtering
#'   and adjusted p-values.
#' 
#' @param d_medians \code{\link[SummarizedExperiment]{SummarizedExperiment}} object 
#'   containing cluster medians (median expression of functional markers), from 
#'   \code{\link{calcMedians}}.
#' 
#' @param d_vals \code{\link[SummarizedExperiment]{SummarizedExperiment}} object 
#'   containing subsetted marker expression values for each cluster-sample combination,
#'   calculated with \code{\link{subsetVals}}.
#' 
#' @param group_IDs Vector or factor of group membership IDs for each sample (e.g. 
#'   diseased vs. healthy, or treated vs. untreated). Vectors are converted to factors 
#'   internally. The first level of the factor will be used as the reference level for 
#'   differential testing. Currently, only two-group comparisons are implemented.
#' 
#' @param paired Whether to perform paired tests. Set to TRUE and provide the 
#'   \code{block_IDs} argument (e.g. patient IDs) to calculate paired tests. Default = 
#'   FALSE.
#' 
#' @param block_IDs Vector or factor of block IDs for samples (e.g. patient ID), required 
#'   for paired tests. Default = NULL.
#' 
#' @param min_cells Filtering parameter. For each cluster-marker combination, samples with
#'   less than \code{min_cells} cells are removed. If \code{paired = TRUE}, both samples 
#'   from a pair are removed if either has less than \code{min_cells} cells. Default = 5.
#' 
#' @param min_samples Filtering parameter. Clusters with less than \code{min_samples} 
#'   samples remaining in either group are removed (for a given marker), with p-value = NA
#'   returned in these cases. Default = 2.
#' 
#' @param n_perm Number of permutations to use for permutation testing (unpaired tests
#'   only). By definition of permutation tests, the minimum p-value is 1/n_perm. Default =
#'   1000 (i.e. minimum p-value = 0.001).
#' 
#' @param n_cores Number of processor cores for parallelized evaluation. Default = all 
#'   available cores.
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
#' @importFrom stats ks.test pchisq
#' @importFrom IHW ihw adj_pvalues
#' @importFrom reshape2 melt
#' @importFrom BiocParallel bplapply MulticoreParam multicoreWorkers
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
#' # subset marker expression values
#' d_vals <- subsetVals(d_se)
#' 
#' 
#' #############################################################################
#' # Test for differential expression (DE) of functional markers within clusters
#' # (method 'diffcyt-KS')
#' #############################################################################
#' 
#' # create block IDs for paired tests (this is a paired data set, so we use 1 block per patient)
#' patient_IDs <- factor(gsub("_(BCRXL|ref)$", "", sample_IDs))
#' patient_IDs <- as.numeric(patient_IDs)
#' patient_IDs
#' 
#' # test for differential expression (DE) of functional markers within clusters
#' # (for demonstration purposes: single core)
#' set.seed(123)
#' res_DE <- testDE_KS(d_counts, d_medians, d_vals, group_IDs, 
#'                     paired = TRUE, block_IDs = patient_IDs)
#' 
#' # show results using 'rowData' accessor function
#' rowData(res_DE)
#' 
#' # sort to show top (most highly significant) cluster-marker combinations first
#' head(rowData(res_DE)[order(rowData(res_DE)$p_adj), ], 10)
#' 
testDE_KS <- function(d_counts, d_medians, d_vals, group_IDs, 
                      paired = FALSE, block_IDs = NULL, 
                      min_cells = 5, min_samples = 2, 
                      n_perm = 1000, n_cores = NULL) {
  
  if (!is.factor(group_IDs)) {
    group_IDs <- factor(group_IDs, levels = unique(group_IDs))
  }
  if (!is.null(block_IDs) & !is.factor(block_IDs)) {
    block_IDs <- factor(block_IDs, levels = unique(block_IDs))
  }
  
  if (paired & is.null(block_IDs)) {
    stop("'block_IDs' argument is required for paired tests")
  }
  
  cluster <- rowData(d_vals)$cluster
  samples <- colData(d_vals)$sample
  func_names <- names(assays(d_vals))
  
  # currently only two-group comparisons possible
  grp_ref <- group_IDs == levels(group_IDs)[1]
  
  
  # ------------------------------------------------------------------------------
  # paired tests: calculate p-values (no permutations or parallelization required)
  # ------------------------------------------------------------------------------
  
  if (paired) {
    p_vals <- matrix(NA, nrow = length(levels(cluster)), ncol = length(func_names))
    rownames(p_vals) <- levels(cluster)
    colnames(p_vals) <- func_names
    
    # clusters 'i' and markers 'j'
    for (j in seq_along(func_names)) {
      for (i in seq_along(cluster)) {
        
        y <- assays(d_vals)[[j]][i, ]
        
        # fix missing values (when zero cells in this cluster-sample combination)
        y <- lapply(y, function(z) {
          if (length(z) == 0) z <- NA
          z
        })
        
        # filtering step 1: minimum number of cells per sample
        n_cells_i <- assays(d_counts)[[1]][i, ]
        keep <- n_cells_i >= min_cells
        if (paired) {
          if (!all(block_IDs[grp_ref] == block_IDs[!grp_ref])) {
            stop(paste0("Block IDs for paired testing are in different order for each group, ", 
                        "which creates problems during filtering. Please rearrange order of samples."))
          }
          keep1 <- grp_ref & keep
          keep2 <- !grp_ref & keep
          keep <- unname(rep(keep1[grp_ref] & keep2[!grp_ref], length(levels(group_IDs))))
        }
        
        y <- y[keep]
        grp_ref_keep <- grp_ref[keep]
        
        # filtering step 2: minimum number of samples per group
        if (length(y) == 0) {
          p_vals[i, j] <- NA
          next
        }
        if ((min(table(grp_ref_keep)) < min_samples) | (length(table(grp_ref_keep)) <= 1)) {
          p_vals[i, j] <- NA
          next
        }
        
        # calculate KS tests (paired)
        y1 <- y[grp_ref_keep]
        y2 <- y[!grp_ref_keep]
        
        if (!(length(y1) == length(y2))) {
          stop("paired tests: number of samples must be the same in both groups")
        }
        
        KS <- rep(NA, length(y1))
        for (z in seq_along(KS)) {
          KS[z] <- ks.test(y1[[z]], y2[[z]], alternative = "two.sided")$p.value
        }
        
        # combine p-values (Fisher's method)
        stat <- -2 * sum(log(KS))
        p_val <- pchisq(stat, df = 2 * length(KS), lower.tail = FALSE)
        
        p_vals[i, j] <- p_val
      }
    }
  }
  
  
  # ----------------------------------
  # unpaired tests: calculate p-values
  # ----------------------------------
  
  if (!paired) {
    
    # function for parallelized evaluation: calculates p-value for cluster 'i' and marker 'j'
    calc_p_val <- function(indices, d_vals, d_counts, min_cells, min_samples, n_perm, grp_ref) {
      
      # extract 'i' and 'j' from 'indices' (data frame)
      i <- unname(unlist(indices)[1])
      j <- unname(unlist(indices)[2])
      
      y <- assays(d_vals)[[j]][i, ]
      
      # fix missing values (when zero cells in this cluster-sample combination)
      y <- lapply(y, function(z) {
        if (length(z) == 0) z <- NA
        z
      })
      
      # filtering step 1: minimum number of cells per sample
      n_cells_i <- assays(d_counts)[[1]][i, ]
      keep <- n_cells_i >= min_cells
      
      y <- y[keep]
      grp_ref_keep <- grp_ref[keep]
      
      # filtering step 2: minimum number of samples per group
      if (length(y) == 0) {
        return(NA)
      }
      if ((min(table(grp_ref_keep)) < min_samples) | (length(table(grp_ref_keep)) <= 1)) {
        return(NA)
      }
      
      # calculate test statistic: maximum KS statistic between any two samples between groups
      y1 <- y[grp_ref_keep]
      y2 <- y[!grp_ref_keep]
      
      KS_stats <- matrix(NA, nrow = length(y1), ncol = length(y2))
      for (z1 in seq_along(y1)) {
        for (z2 in seq_along(y2)) {
          KS_stats[z1, z2] <- ks.test(y1[[z1]], y2[[z2]], alternative = "two.sided")$statistic
        }
      }
      KS_min <- min(KS_stats)
      
      # permutation null distributions: permute group labels
      KS_perm <- rep(NA, n_perm)
      
      for (b in seq_len(n_perm)) {
        perm <- sample(grp_ref_keep)
        y1 <- y[perm]
        y2 <- y[!perm]
        
        KS_stats <- matrix(NA, nrow = length(y1), ncol = length(y2))
        for (z1 in seq_along(y1)) {
          for (z2 in seq_along(y2)) {
            KS_stats[z1, z2] <- ks.test(y1[[z1]], y2[[z2]], alternative = "two.sided")$statistic
          }
        }
        KS_perm[b] <- min(KS_stats)
      }
      
      # final p-value
      p_val <- mean(KS_min < KS_perm)
      p_val
    }
    
    
    # BiocParallel parameters
    if (is.null(n_cores)) n_cores <- multicoreWorkers()
    bpparam <- MulticoreParam(workers = n_cores)
    
    # set up grid of indices (i, j) for parallelization
    indices <- expand.grid(seq_along(cluster), seq_along(func_names))
    indices <- split(indices, seq_len(nrow(indices)))
    
    # calculate p-values (parallelized across clusters 'i' and markers 'j')
    p_vals <- bplapply(indices, calc_p_val, d_vals, d_counts, min_cells, min_samples, 
                       n_perm, grp_ref, BPPARAM = bpparam)
    
    p_vals <- matrix(unlist(p_vals), nrow = length(levels(cluster)), ncol = length(func_names))
    rownames(p_vals) <- levels(cluster)
    colnames(p_vals) <- func_names
  }
  
  
  # --------------
  # return results
  # --------------
  
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
    lapply(as.list(assays(d_medians)[func_names]), function(a) a[cluster, ])
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


