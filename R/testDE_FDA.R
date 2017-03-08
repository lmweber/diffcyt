#' Test for differential expression of functional markers (method: "diffcyt-FDA")
#' 
#' Calculate tests for differential expression of functional markers within clusters 
#' (method: "diffcyt-FDA")
#' 
#' The "diffcyt-FDA" methodology uses techniques from functional data analysis (FDA) to 
#' model the functional marker expression signals within each cluster. The differential 
#' expression tests compare the FDA-modeled expression profiles within clusters between 
#' samples in the two groups (e.g. diseased vs. healthy).
#' 
#' Permutation tests are used to calculate p-values.
#' 
#' By default, the tests are weighted by the number of cells per cluster-sample
#' combination, representing the relative uncertainty in calculating each curve (i.e. the
#' ECDF for each cluster-sample combination). Note that these are relative weights within
#' each cluster only. Weights can also be disabled for (much) faster runtime
#' (\code{weights = FALSE}), but then the results will not account for this source of
#' uncertainty.
#' 
#' Two levels of filtering are included by default. First, for each cluster-marker
#' combination, samples with less than \code{min_cells} cells are removed. Second, any
#' clusters with less than \code{min_samples} samples in either condition are removed (for
#' this marker); p-value = NA is returned in these cases.
#' 
#' We use the \code{\link[fda]{fda}} package (Ramsay et al. 2014) for the FDA modeling 
#' steps; with modifications to the permutation t-testing function to allow weights and 
#' improve runtime.
#' 
#' The \code{\link[BiocParallel]{BiocParallel}} package is used for parallelized 
#' evaluation on multi-core systems, to speed up runtime.
#' 
#' Alternative methodologies for testing for differential expression of functional markers
#' are available in the functions \code{testDE_med} and \code{testDE-KS}.
#' 
#' 
#' @param d_ecdfs \code{\link[SummarizedExperiment]{SummarizedExperiment}} object 
#'   containing empirical cumulative distribution functions (ECDFs), calculated with 
#'   \code{\link{calcECDFs}}.
#' 
#' @param d_clus \code{\link[SummarizedExperiment]{SummarizedExperiment}} object 
#'   containing cluster medians (median expression of functional markers) and cluster 
#'   frequencies (number of cells), from \code{\link{calcMediansAndFreq}}. The cluster 
#'   frequencies are used as weights for the permutation t-tests (the cluster medians are 
#'   not used here).
#' 
#' @param group Factor containing group membership for each sample (for example, diseased 
#'   vs. healthy), for differential comparisons and statistical tests.
#' 
#' @param weighted Whether to include weights (per cluster-sample combination) for 
#'   weighted permutation tests (see details below). Without weights, runtime is much 
#'   faster, but result do not account for uncertainty due to different numbers of cells 
#'   per sample (within each cluster). Default = TRUE.
#' 
#' @param min_cells Filtering parameter. For each cluster-marker combination, samples with
#'   less than \code{min_cells} cells are removed. Default = 6.
#' 
#' @param min_samples Filtering parameter. Clusters with less than \code{min_samples} 
#'   samples in either group (condition) are removed (for a given marker); p-value = NA in
#'   these cases. Default = 2.
#' 
#' @param n_perm Number of permutations to use for permutation testing. Default = 5000 
#'   (i.e. minimum possible p-value = 0.0002).
#' 
#' @param n_cores Number of processor cores for parallelized evaluation. Default = all
#'   available cores.
#' 
#' 
#' @return Returns a data frame containing cluster labels, functional marker names, and 
#'   p-values summarizing the evidence for differential expression. Rows are sorted by 
#'   p-value (smallest first).
#' 
#' 
#' @importFrom SummarizedExperiment assays rowData colData
#' @importFrom fda Data2fd
#' @importFrom reshape2 melt
#' @importFrom BiocParallel bplapply MulticoreParam multicoreWorkers
#' 
#' @export
#' 
#' @examples
#' library(flowCore)
#' library(limma)
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
#' # indices of all marker columns, lineage markers, and functional markers
#' # (see Table 1 in Bruggner et al. 2014)
#' marker_cols <- c(3:4, 7:9, 11:19, 21:22, 24:26, 28:31, 33)
#' lineage_cols <- c(3:4, 9, 11,12,14, 21, 29, 31, 33)
#' functional_cols <- setdiff(marker_cols, lineage_cols)
#' 
#' # prepare data
#' d_se <- prepareData(d_input, sample_IDs, group_IDs, marker_cols, lineage_cols, functional_cols)
#' 
#' # transform data
#' d_se <- transformData(d_se, cofactor = 5)
#' 
#' # generate clusters (note: using small 3x3 SOM grid for demonstration purposes)
#' d_se <- generateClusters(d_se, cols_to_use = lineage_cols, xdim = 3, ydim = 3, 
#'                          seed = 123, plot = FALSE)
#' 
#' # calculate cluster medians and frequencies
#' d_clus <- calcMediansAndFreq(d_se)
#' 
#' # calculate ECDFs
#' d_ecdfs <- calcECDFs(d_se)
#' 
#' 
#' #########################################################################################
#' # (3) Test for differential expression (DE) of functional markers: method "diffcyt-FDA" #
#' #########################################################################################
#' 
#' # re-level factor to use "ref" as base level
#' group <- factor(group_IDs, levels = c("ref", "BCRXL"))
#' 
#' # note: using no weights, small number of permutations, and single core for 
#' # demonstration purposes
#' res_DE_FDA <- testDE_FDA(d_ecdfs, d_clus, group, weighted = FALSE, 
#'                          n_perm = 100, n_cores = 1)
#' 
#' # (note: this is a small example data set used for demonstration purposes only; results
#' # are not biologically meaningful)
#' head(res_DE_FDA)
#' 
testDE_FDA <- function(d_ecdfs, d_clus, group, weighted = TRUE, 
                       min_cells = 6, min_samples = 2, 
                       n_perm = 5000, n_cores = NULL) {
  
  if (!is.factor(group)) group <- factor(group, levels = unique(group))
  
  # note: assumes same resolution for all ECDFs
  resolution <- length(assays(d_ecdfs)[[1]][1, 1][[1]])
  
  clus <- rowData(d_ecdfs)$cluster
  smp <- colData(d_ecdfs)$sample
  func_markers <- names(assays(d_ecdfs))
  
  # number of cells per cluster-sample combination: to use as weights
  # [to do: include check that column order matches groups]
  if (weighted) {
    weights <- assays(d_clus)[["n_cells"]]
  } else {
    weights <- NULL
  }
  
  # [to do: generalize for different experimental designs]
  grp <- group == levels(group)[1]
  
  argvals <- seq_len(resolution)
  
  # function for parallelized evaluation: calculates p-value for cluster 'i' and marker 'j'
  calc_p_val <- function(indices, d_ecdfs, d_clus, min_cells, min_samples, resolution, 
                         smp, argvals, grp, weighted, weights, n_perm) {
    
    # extract 'i' and 'j' from 'indices' (data frame)
    i <- unname(unlist(indices)[1])
    j <- unname(unlist(indices)[2])
    
    y <- assays(d_ecdfs)[[j]][i, ]
    
    # fix missing values (when zero cells in this cluster-sample combination)
    # [to do: move this to 'calcECDFs()']
    y <- lapply(y, function(yy) {
      if (length(yy) == 0) yy <- rep(0, resolution)
      yy
    })
    
    # filtering step 1: minimum number of cells per sample
    n_cells_i <- assays(d_clus)[["n_cells"]][i, ]
    keep <- n_cells_i >= min_cells
    y <- y[keep]
    grp <- grp[keep]
    
    if (weighted) {
      weights1 <- weights[, grp]
      weights2 <- weights[, !grp]
    }
    
    # filtering step 2: minimum number of samples per condition
    if ((min(table(grp)) <= min_samples) | (length(table(grp)) <= 1)) {
      return(NA)
    }
    
    # prepare 'fda' objects
    y <- matrix(unlist(y), ncol = length(grp))
    fd1 <- Data2fd(argvals = argvals, y = y[, grp])
    fd2 <- Data2fd(argvals = argvals, y = y[, !grp])
    
    # note: keeping p-values only (discarding all other results)
    if (weighted) {
      p_val <- .tperm.fd_wtd_fast(fd1, fd2, weights1[i, ], weights2[i, ], nperm = n_perm, plotres = FALSE)$pval
    } else {
      p_val <- .tperm.fd_fast(fd1, fd2, nperm = n_perm, plotres = FALSE)$pval
    }
    
    p_val
  }
  
  
  # BiocParallel parameters
  if (is.null(n_cores)) n_cores <- multicoreWorkers()
  bpparam <- MulticoreParam(workers = n_cores)
  
  # grid of indices for parallelization
  indices <- expand.grid(seq_along(clus), seq_along(func_markers))
  indices <- split(indices, seq_len(nrow(indices)))
  
  # calculate p-values (parallelized across clusters 'i' and markers 'j')
  p_vals <- bplapply(indices, calc_p_val, d_ecdfs, d_clus, min_cells, min_samples, 
                     resolution, smp, argvals, grp, weighted, weights, n_perm, 
                     BPPARAM = bpparam)
  
  p_vals <- matrix(unlist(p_vals), nrow = length(clus), ncol = length(func_markers))
  rownames(p_vals) <- clus
  colnames(p_vals) <- func_markers
  
  # sort to show lowest p-values (across functional markers and clusters) at the top
  res <- as.data.frame(p_vals)
  res$cluster <- rownames(res)
  res <- melt(res, id.vars = "cluster", variable.name = "marker", value.name = "p_val")
  
  res <- res[order(res$p_val), ]
  
  res
}


