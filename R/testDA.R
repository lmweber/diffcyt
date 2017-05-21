#' Test for differential abundance
#' 
#' Calculate tests for differential abundance of clusters
#' 
#' Calculates tests for differential abundance (differential frequencies) of clusters, 
#' using empirical Bayes moderation of cluster variances to improve power.
#' 
#' We use the \code{\link[limma]{limma}} package (Ritchie et al. 2015, \emph{Nucleic Acids
#' Research}) to calculate the empirical Bayes moderated tests. Empirical Bayes methods 
#' improve statistical power by sharing information on variability (i.e. variance across 
#' samples for a single cluster) between clusters.
#' 
#' Since count data are often heteroscedastic, we use the  \code{\link[limma]{voom}} 
#' method (Law et al. 2014, \emph{Genomie Biology}) to transform the raw cluster cell 
#' counts and estimate observation-level weights to stabilize the mean-variance 
#' relationship. Diagnostic plots are shown if \code{plot = TRUE}.
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
#' For paired data (e.g. treated vs. untreated samples from the same patient), paired 
#' tests can be performed, which improves the statistical power of the differential tests.
#' Paired tests are performed by setting \code{paired = TRUE} and providing the 
#' \code{block_IDs} argument (e.g. one block per patient).
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
#' @param paired Whether to perform paired tests. Set to TRUE and provide the 
#'   \code{block_IDs} argument (e.g. patient IDs) to calculate paired tests. Default = 
#'   FALSE.
#' 
#' @param block_IDs Vector or factor of block IDs for samples (e.g. patient ID), required 
#'   for paired tests. Default = NULL.
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
#' ##################################################
#' # Test for differential abundance (DA) of clusters
#' ##################################################
#' 
#' # create block IDs for paired tests (this is a paired data set, so we use 1 block per patient)
#' patient_IDs <- factor(gsub("_(BCRXL|ref)$", "", sample_IDs))
#' patient_IDs <- as.numeric(patient_IDs)
#' patient_IDs
#' 
#' # test for differential abundance (DA) of clusters
#' res_DA <- testDA(d_counts, group_IDs, paired = TRUE, block_IDs = patient_IDs)
#' 
#' # show results using 'rowData' accessor function
#' rowData(res_DA)
#' 
testDA <- function(d_counts, group_IDs, paired = FALSE, block_IDs = NULL, plot = FALSE, path = ".") {
  
  if (!is.factor(group_IDs)) group_IDs <- factor(group_IDs, levels = unique(group_IDs))
  
  if (paired & is.null(block_IDs)) {
    stop("'block_IDs' argument is required for paired tests")
  }
  
  counts <- assays(d_counts)[[1]]
  cluster <- rowData(d_counts)$cluster
  
  # model matrix
  mm <- model.matrix(~ group_IDs)
  
  # voom transformation and weights
  if (plot) {
    pdf(file.path(path, "mean_variance_pre_voom.pdf"), width = 6, height = 6)
    v <- voom(counts, design = mm, plot = TRUE)
    dev.off()
  } else {
    v <- voom(counts, design = mm, plot = FALSE)
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
  top <- topTable(efit, coef = 2, number = Inf, adjust.method = "BH", sort.by = "none")
  if (!all(rownames(top) == cluster)) {
    stop("cluster labels do not match")
  }
  
  res_DA <- d_counts
  
  rowData(res_DA) <- cbind(rowData(res_DA), top)
  
  # also store additional sample information in 'colData'
  colData(res_DA) <- cbind(colData(res_DA), data.frame(group_IDs), data.frame(block_IDs))
  
  res_DA
}


