#' Test for differential expression within clusters (method: 'diffcyt-med')
#' 
#' Calculate tests for differential expression of functional markers within clusters 
#' (method: 'diffcyt-med'), using empirical Bayes moderation of cluster variances to 
#' improve power.
#' 
#' The 'diffcyt-med' methodology uses median marker expression to characterize the signal
#' of interest within each cluster. The differential expression tests compare the median 
#' marker expression within clusters between samples in the two groups (e.g. diseased vs.
#' healthy).
#' 
#' To stabilize the mean-variance relationship, we use the \code{\link[limma]{voom}} 
#' method (Law et al. 2014, \emph{Genomie Biology}) to transform the expression values and
#' estimate observation-level weights. Diagnostic plots are shown if \code{plot = TRUE}.
#' 
#' Filtering: Clusters are kept for testing if there are at least \code{min_cells} cells
#' per sample in at least \code{min_samples} samples in either condition. Filtered
#' clusters are removed from differential expression testing for all markers.
#' 
#' Empirical Bayes methods are used to share information on variability (i.e. variance 
#' across samples for a single cluster) between clusters, to improve statistical power. We
#' use the \code{\link[limma]{limma}} package (Ritchie et al. 2015, \emph{Nucleic Acids 
#' Research}) to calculate the empirical Bayes moderated tests.
#' 
#' Alternative methodologies for testing for differential expression within clusters are 
#' available in the functions \code{\link{testDE_FDA}}, \code{\link{testDE_KS}}, and 
#' \code{\link{testDE_LM}}.
#' 
#' 
#' @param d_medians \code{\link[SummarizedExperiment]{SummarizedExperiment}} object 
#'   containing cluster medians (median expression of functional markers), from 
#'   \code{\link{calcMedians}}.
#' 
#' @param d_counts \code{\link[SummarizedExperiment]{SummarizedExperiment}} object 
#'   containing cluster cell counts, from \code{\link{calcCounts}}.
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
#' @param plot Whether to save 'voom' diagnostic plots. Default = FALSE.
#' 
#' @param path Path to save diagnostic plots, if \code{plot = TRUE}. Default = current 
#'   working directory.
#' 
#' 
#' @return Returns new \code{\link[SummarizedExperiment]{SummarizedExperiment}} object, 
#'   where rows = clusters, columns = samples. Rows (clusters) are repeated for each 
#'   functional marker (i.e. the sheets or 'assays' from the previous \code{d_medians}
#'   object are stacked into a single matrix). Differential test results are stored in the
#'   'rowData' slot. Results include p-values and adjusted p-values from the
#'   \code{\link[limma]{limma}} empirical Bayes moderated tests, which can be used to rank
#'   clusters (across all functional markers) by evidence for differential expression. The
#'   results can be accessed with the \code{\link[SummarizedExperiment]{rowData}} accessor
#'   function.
#' 
#' 
#' @importFrom SummarizedExperiment assays rowData colData
#' @importFrom limma voom duplicateCorrelation lmFit eBayes plotSA topTable
#' @importFrom stats model.matrix
#' @importFrom methods as is
#' @importFrom grDevices pdf
#' @importFrom graphics plot
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
#' # set group reference level for differential testing
#' group_IDs <- factor(group_IDs, levels = c("ref", "BCRXL"))
#' group_IDs
#' 
#' # indices of all marker columns, lineage markers, and functional markers
#' # (see Table 1 in Bruggner et al. 2014)
#' cols_markers <- c(3:4, 7:9, 11:19, 21:22, 24:26, 28:31, 33)
#' cols_lineage <- c(3:4, 9, 11,12,14, 21, 29, 31, 33)
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
#' 
#' #############################################################################
#' # Test for differential expression (DE) of functional markers within clusters
#' # (method 'diffcyt-med')
#' #############################################################################
#' 
#' # create block IDs for paired tests (this is a paired data set, so we use 1 block per patient)
#' patient_IDs <- factor(gsub("_(BCRXL|ref)$", "", sample_IDs))
#' patient_IDs <- as.numeric(patient_IDs)
#' patient_IDs
#' 
#' # test for differential expression (DE) of functional markers within clusters
#' res_DE <- testDE_med(d_medians, d_counts, group_IDs, paired = TRUE, block_IDs = patient_IDs)
#' 
#' # show results using 'rowData' accessor function
#' rowData(res_DE)
#' 
testDE_med <- function(d_medians, d_counts, group_IDs, paired = FALSE, block_IDs = NULL, 
                       min_cells = 5, min_samples = NULL, plot = FALSE, path = ".") {
  
  if (!is.factor(group_IDs)) group_IDs <- factor(group_IDs, levels = unique(group_IDs))
  
  if (paired & is.null(block_IDs)) {
    stop("'block_IDs' argument is required for paired tests")
  }
  
  if (is.null(min_samples)) {
    min_samples <- min(table(group_IDs)) - 1
  }
  
  counts <- assays(d_counts)[[1]]
  cluster <- rowData(d_counts)$cluster
  
  # filtering
  grp <- group_IDs == levels(group_IDs)[1]
  tf <- counts >= min_cells
  ix_keep <- (rowSums(tf[, grp]) >= min_samples) & (rowSums(tf[, !grp]) >= min_samples)
  
  counts <- counts[ix_keep, ]
  cluster <- cluster[ix_keep]
  
  # remove any remaining rows with zeros (required by voom)
  ix_zeros <- apply(counts_tmp, 1, function(r) any(r == 0))
  
  counts <- counts[!ix_zeros, ]
  cluster <- cluster[!ix_zeros]
  
  ix_meds <- cluster
  
  # extract medians and create concatenated matrix
  func_names <- names(assays(d_medians))
  meds <- do.call("rbind", {
    lapply(as.list(assays(d_medians)[func_names]), function(a) a[ix_meds, ])
  })
  
  # functional marker names and cluster labels (for SummarizedExperiment rowData)
  clus <- rep(cluster, length(func_names))
  func <- rep(func_names, each = length(cluster))
  
  if (!(length(clus) == nrow(meds))) {
    stop("Check filtering: length of cluster labels vector does not match medians matrix")
  }
  if (!(length(func) == nrow(meds))) {
    stop("Check filtering: length of functional marker names vector does not match medians matrix")
  }
  
  # model matrix
  mm <- model.matrix(~ group_IDs)
  
  # voom transformation and weights
  if (plot) {
    pdf(file.path(path, "mean_variance_pre_voom.pdf"), width = 6, height = 6)
    v <- voom(meds, design = mm, plot = TRUE)
    dev.off()
  } else {
    v <- voom(meds, design = mm, plot = FALSE)
  }
  
  # fit linear models
  if (paired) {
    dupcor <- duplicateCorrelation(v, design = mm, block = block_IDs)
    vfit <- lmFit(v, design = mm, block = block_IDs, correlation = dupcor$consensus.correlation)
  } else {
    vfit <- lmFit(v, design = mm)
  }
  
  # calculate empirical Bayes moderated tests
  efit <- eBayes(vfit)
  
  if (plot) {
    pdf(file.path(path, "mean_variance_post_voom.pdf"), width = 6, height = 6)
    plotSA(efit)
    dev.off()
  }
  
  # return new 'SummarizedExperiment' object with results stored in 'rowData'
  col_data <- cbind(colData(d_medians), data.frame(group_IDs), data.frame(block_IDs))
  
  top <- topTable(efit, coef = 2, number = Inf, adjust.method = "BH", sort.by = "none")
  if (!all(top$ID %in% cluster) | !(nrow(top) == nrow(meds))) {
    stop("cluster labels do not match")
  }
  
  row_data <- cbind(data.frame(cluster = clus, marker = func), top)
  
  res_DE <- SummarizedExperiment(meds, rowData = row_data, colData = col_data)
  
  res_DE
}


