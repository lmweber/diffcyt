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
#' We use the \code{\link[fda]{fda}} package (Ramsay et al. 2014) for the FDA modeling 
#' steps.
#' 
#' Alternative methodologies for testing for differential expression of functional markers
#' are available in the functions \code{testDE_med} and \code{testDE-KS}.
#' 
#' 
#' @param d_ecdfs \code{\link[SummarizedExperiment]{SummarizedExperiment}} object 
#'   containing empirical cumulative distribution functions (ECDFs), calculated with 
#'   \code{\link{calcECDFs}}.
#' 
#' @param group Factor containing group membership for each sample (for example, diseased 
#'   vs. healthy), for differential comparisons and statistical tests.
#' 
#' @param n_perm Number of permutations to use for permutation testing. Default = 5000 
#'   (i.e. minimum possible p-value = 0.0002).
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
#' # generate clusters (small 10x10 SOM grid due to small size of example data set)
#' d_se <- generateClusters(d_se, cols_to_use = lineage_cols, xdim = 10, ydim = 10, 
#'                          seed = 123, plot = FALSE)
#' 
#' # calculate cluster medians and frequencies
#' d_clus <- calcMediansAndFreq(d_se)
#' 
#' # calculate ECDFs
#' d_ecdfs <- calcECDFs(d_se)
#' 
#' # (1) test for differentially abundant (DA) clusters
#' group <- factor(group_IDs, levels = c("ref", "BCRXL"))  # re-level factor to use "ref" as base level
#' res_DA <- testDA(d_clus, group)
#' topTable(res_DA, number = 6)
#' 
#' # plot top DA clusters
#' # note there is no evidence for DA in this example data set (data set is too small)
#' # plotTopDAClusters(res_DA)
#' 
#' # (3) test for differential expression of functional markers: method "diffcyt-FDA"
#' # (note: using small number of permutations for demonstration purposes)
#' res_DE_FDA <- testDE_FDA(d_ecdfs, group, n_perm = 100)
#' head(res_DE_FDA)
testDE_FDA <- function(d_ecdfs, group, n_perm = 5000) {
  
  if (!is.factor(group)) group <- factor(group, levels = unique(group))
  
  # note: assumes same resolution for all ECDFs
  resolution <- length(assays(d_ecdfs)[[1]][1, 1][[1]])
  
  clus <- rowData(d_ecdfs)$cluster
  smp <- colData(d_ecdfs)$sample
  func_markers <- names(assays(d_ecdfs))
  
  # set up matrix for p-values
  p_vals <- matrix(NA, nrow = length(clus), ncol = length(func_markers))
  rownames(p_vals) <- clus
  colnames(p_vals) <- func_markers
  
  # [to do: generalize for different experimental designs]
  grp <- group == levels(group)[1]
  
  argvals <- 1:resolution
  
  # calculate p-values
  # [to do: include some filtering: no. of cells, no. of samples]
  for (j in seq_along(func_markers)) {
    
    assays_j <- assays(d_ecdfs)[[j]]
    
    for (i in seq_along(clus)) {
      cat(".")
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
        .tperm.fd_fast(fd1, fd2, nperm = n_perm, plotres = FALSE)$pval
      }, error = function(e) e)
      
      if (is(p_val, "simpleError")) p_val <- NA
      
      p_vals[i, j] <- p_val
    }
  }
  
  # sort to show lowest p-values (across all cluster-sample combinations) at the top
  res <- as.data.frame(p_vals)
  res$cluster <- rownames(res)
  res <- melt(res, id.vars = "cluster", variable.name = "marker", value.name = "p_val")
  
  res <- res[order(res$p_val), ]
  
  res
}


