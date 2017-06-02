#' Test for differential abundance
#' 
#' Calculate tests for differential abundance of clusters
#' 
#' Calculates tests for differential abundance (differential cell frequencies) of
#' clusters, using empirical Bayes moderation of cluster variances to improve power.
#' 
#' We use the \code{\link[limma]{limma}} package (Ritchie et al. 2015, \emph{Nucleic Acids
#' Research}) to calculate the empirical Bayes moderated tests. Empirical Bayes methods 
#' improve statistical power by sharing information on variability (i.e. variance across 
#' samples for a single cluster) between clusters.
#' 
#' Since count data are often heteroscedastic, we use the  \code{\link[limma]{voom}} 
#' method (Law et al. 2014, \emph{Genome Biology}) to transform the raw cluster cell 
#' counts and estimate observation-level weights to stabilize the mean-variance 
#' relationship. Diagnostic plots are shown if \code{plot = TRUE}.
#' 
#' Filtering: Clusters are kept for testing if there are at least \code{min_cells} cells 
#' per sample in at least \code{min_samples} samples in at least one condition. The
#' \code{voom} diagnostic plots can be used to optimize the level of filtering: the
#' mean-variance trend prior to transformation (sqrt standard deviation vs. log2 count
#' size) should be strictly monotonically decreasing; if there is an initial increasing
#' trend, the amount of filtering should be increased (Law et al. 2014, \emph{Genome
#' Biology}).
#' 
#' The group membership IDs for each sample (e.g. diseased vs. healthy, or treated vs. 
#' untreated) are specified with the \code{group_IDs} argument. The group IDs should be 
#' checked carefully to ensure they are in the same order as the samples (columns) in the 
#' data object \code{d_counts}.
#' 
#' The \code{group_IDs} argument can be provided as a vector or factor. Vectors are 
#' converted to factors internally. For two-group comparisons, the first level of the 
#' factor will be used as the reference level for the differential tests. To specify a 
#' different reference level, use \code{relevel} to re-order the levels (so that the 
#' reference level is first).
#' 
#' Currently, only two-group comparisons are possible. More complex comparisons 
#' (contrasts) will be implemented in a future version.
#' 
#' 
#' @param d_counts \code{\link[SummarizedExperiment]{SummarizedExperiment}} object 
#'   containing cluster cell counts, from \code{\link{calcCounts}}.
#' 
#' @param group_IDs Vector or factor of group membership IDs for each sample (e.g. 
#'   diseased vs. healthy, or treated vs. untreated). Vectors are converted to factors 
#'   internally. The first level of the factor will be used as the reference level for 
#'   differential testing. Currently, only two-group comparisons are implemented.
#' 
#' @param min_cells Filtering parameter. Default = 5. Clusters are kept if there are at 
#'   least \code{min_cells} cells per sample in at least \code{min_samples} samples in at
#'   least one condition.
#' 
#' @param min_samples Filtering parameter. Default = \code{n - 1}, where \code{n} = number
#'   of replicates in smallest group. Clusters are kept if there are at least 
#'   \code{min_cells} cells per sample in at least \code{min_samples} samples in at least
#'   one condition.
#' 
#' @param plot Whether to save 'voom' diagnostic plots. Default = FALSE.
#' 
#' @param path Path to save diagnostic plots, if \code{plot = TRUE}. Default = current 
#'   working directory.
#' 
#' 
#' @return Returns new \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'   with differential test results for each cluster stored in the 'rowData' slot. Results
#'   include p-values and adjusted p-values from the \code{\link[limma]{limma}} empirical
#'   Bayes moderated tests, which can be used to rank clusters by evidence for
#'   differential abundance. The results can be accessed with the 
#'   \code{\link[SummarizedExperiment]{rowData}} accessor function.
#' 
#' 
#' @importFrom SummarizedExperiment assays rowData 'rowData<-' colData 'colData<-'
#' @importFrom limma voom lmFit eBayes plotSA topTable
#' @importFrom stats model.matrix
#' @importFrom methods as is
#' @importFrom grDevices pdf
#' @importFrom graphics plot
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
#' ################################################
#' # Test for differentially abundant (DA) clusters
#' ################################################
#' 
#' # test for differentially abundant (DA) clusters
#' res_DA <- testDA(d_counts, group_IDs)
#' 
#' # show results using 'rowData' accessor function
#' rowData(res_DA)
#' 
#' # sort to show top (most highly significant) clusters first
#' head(rowData(res_DA)[order(rowData(res_DA)$adj.P.Val), ], 10)
#' 
testDA <- function(d_counts, group_IDs, 
                   min_cells = 5, min_samples = NULL, 
                   plot = FALSE, path = ".") {
  
  if (!is.factor(group_IDs)) group_IDs <- factor(group_IDs, levels = unique(group_IDs))
  
  if (is.null(min_samples)) {
    min_samples <- min(table(group_IDs)) - 1
  }
  
  counts <- assays(d_counts)[[1]]
  cluster <- rowData(d_counts)$cluster
  
  # filtering
  grp <- group_IDs == levels(group_IDs)[1]
  tf <- counts >= min_cells
  ix_keep <- (rowSums(tf[, grp]) >= min_samples) | (rowSums(tf[, !grp]) >= min_samples)
  
  counts <- counts[ix_keep, ]
  cluster <- cluster[ix_keep]
  
  # model matrix
  mm <- model.matrix(~ group_IDs)
  
  # voom transformation and weights
  if (plot) {
    pdf(file.path(path, "testDA_mean_var_pre_voom.pdf"), width = 6, height = 6)
    v <- voom(counts, design = mm, plot = TRUE)
    dev.off()
  } else {
    v <- voom(counts, design = mm, plot = FALSE)
  }
  
  # fit linear models
  vfit <- lmFit(v, design = mm)
  
  # calculate empirical Bayes moderated tests
  efit <- eBayes(vfit)
  
  if (plot) {
    pdf(file.path(path, "testDA_mean_var_post_voom.pdf"), width = 6, height = 6)
    plotSA(efit)
    dev.off()
  }
  
  # return new 'SummarizedExperiment' object with results stored in 'rowData'
  
  top <- topTable(efit, coef = 2, number = Inf, adjust.method = "BH", sort.by = "none")
  if (!all(rownames(top) == cluster)) {
    stop("cluster labels do not match")
  }
  
  # fill in missing rows (filtered clusters) with NAs
  row_data <- as.data.frame(matrix(NA, nrow = nlevels(cluster), ncol = ncol(top)))
  colnames(row_data) <- colnames(top)
  cluster_nm <- as.numeric(cluster)
  row_data[cluster_nm, ] <- top
  
  row_data <- cbind(data.frame(cluster = as.numeric(levels(cluster))), row_data)
  
  # also store additional sample information in 'colData'
  col_data <- cbind(colData(d_counts), data.frame(group_IDs))
  
  res_DA <- d_counts
  
  rowData(res_DA) <- row_data
  colData(res_DA) <- col_data
  
  res_DA
}


