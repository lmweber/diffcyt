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
#' The tests are weighted by the number of cells per cluster-sample combination,
#' representing the relative uncertainty in calculating each curve (i.e. the ECDF for each
#' cluster-sample combination). Note that these are relative weights within each cluster
#' only.
#' 
#' We use the \code{\link[fda]{fda}} package (Ramsay et al. 2014) for the FDA modeling 
#' steps; with modifications to the permutation t-testing function to allow weighting and
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
#' # note: using small number of permutations for demonstration purposes
#' res_DE_FDA <- testDE_FDA(d_ecdfs, d_clus, group, n_perm = 10)
#' 
#' # (note: this is a small example data set used for demonstration purposes only; results
#' # are not biologically meaningful)
#' head(res_DE_FDA)
#' 
testDE_FDA <- function(d_ecdfs, d_clus, group, n_perm = 5000, n_cores = NULL) {
  
  if (!is.factor(group)) group <- factor(group, levels = unique(group))
  
  # note: assumes same resolution for all ECDFs
  resolution <- length(assays(d_ecdfs)[[1]][1, 1][[1]])
  
  clus <- rowData(d_ecdfs)$cluster
  smp <- colData(d_ecdfs)$sample
  func_markers <- names(assays(d_ecdfs))
  
  # number of cells per cluster-sample combination: to use as weights
  # [to do: include check that column order matches groups]
  weights <- assays(d_clus)[["n_cells"]]
  
  # set up matrix for p-values
  p_vals <- matrix(NA, nrow = length(clus), ncol = length(func_markers))
  rownames(p_vals) <- clus
  colnames(p_vals) <- func_markers
  
  # [to do: generalize for different experimental designs]
  grp <- group == levels(group)[1]
  
  argvals <- 1:resolution
  
  weights1 <- weights[, grp]
  weights2 <- weights[, !grp]
  
  
  # function for parallelized evaluation: calculates p-value for cluster 'i'
  eval_iter_i <- function(i, assays_j, resolution, smp, argvals, grp, weights1, weights2, n_perm) {
    
    y <- assays_j[i, ]
    
    # fix missing values (when zero cells in this cluster-sample combination)
    # [to do: move this to 'calcECDFs()']
    y <- lapply(y, function(yy) {
      if (length(yy) == 0) yy <- rep(0, resolution)
      yy
    })
    
    y <- matrix(unlist(y), ncol = length(smp))
    
    fd1 <- Data2fd(argvals = argvals, y = y[, grp])
    fd2 <- Data2fd(argvals = argvals, y = y[, !grp])
    
    # note: keeping p-values only (discarding all other results)
    # note: error handling: may return errors for some iterations; return NA in these cases
    p_val <- tryCatch({
      .tperm.fd_wtd_fast(fd1, fd2, weights1[i, ], weights2[i, ], nperm = n_perm, plotres = FALSE)$pval
    }, error = function(e) e)
    
    if (is(p_val, "simpleError")) p_val <- NA
    
    p_val
  }
  
  
  # BiocParallel parameters
  if (is.null(n_cores)) n_cores <- multicoreWorkers()
  bpparam <- MulticoreParam(workers = n_cores)
  
  # calculate p-values (parallelized for each marker 'j')
  # [to do: include some filtering: no. of cells, no. of samples]
  for (j in seq_along(func_markers)) {
    
    assays_j <- assays(d_ecdfs)[[j]]
    
    p_vals[, j] <- unlist(
      bplapply(seq_along(clus), eval_iter_i, assays_j, resolution, smp, argvals, grp, weights1, weights2, n_perm, 
               BPPARAM = bpparam)
    )
    
    cat("marker", j, "complete\n")
  }
  
  
  # sort to show lowest p-values (across functional markers and clusters) at the top
  res <- as.data.frame(p_vals)
  res$cluster <- rownames(res)
  res <- melt(res, id.vars = "cluster", variable.name = "marker", value.name = "p_val")
  
  res <- res[order(res$p_val), ]
  
  res
}


