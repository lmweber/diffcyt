#' Test for differential expression of functional markers (method: "medians")
#' 
#' Calculate tests for differential expression of functional markers within clusters 
#' (method: "medians"), using empirical Bayes moderation of variances to improve power.
#' 
#' The "diffcyt-med" methodology uses median functional marker expression to characterize
#' the signal within each cluster. The differential expression tests compare the median
#' functional marker expression within clusters between samples in the two groups (e.g.
#' diseased vs. healthy).
#' 
#' The number of cells per cluster per sample is used to weight the tests, representing
#' the uncertainty in calculating each median value.
#' 
#' Empirical Bayes methods are used to share information on variability across clusters, 
#' improving power. We use the \code{\link[limma]{limma}} package (Ritchie et al. 2015, 
#' \emph{Nucleic Acids Research}) to calculate the empirical Bayes moderated tests.
#' 
#' Alternative methodologies for testing for differential expression of functional markers
#' are available in the functions \code{testDE-FDA} and \code{testDE-KS}.
#' 
#' 
#' @param d_clus \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'   containing cluster medians (median expression of functional markers) and cluster
#'   frequencies (number of cells), from \code{\link{calcMediansAndFreq}}.
#' 
#' @param group Factor containing group membership for each sample (for example, diseased
#'   vs. healthy), for differential comparisons and statistical tests.
#' 
#' @param plot Whether to save plot. Default = TRUE.
#' 
#' @param path Path to save plot.
#' 
#' @param filename Filename for plot.
#' 
#' 
#' @return Returns fitted model objects from \code{\link[limma]{limma}} for each cluster, 
#'   including p-values from the empirical Bayes moderated tests for differential 
#'   expression of functional markers. The p-values can be accessed with the function 
#'   \code{\link[limma]{topTable}} from the \code{limma} package.
#' 
#' 
#' @importFrom limma lmFit eBayes
#' @importFrom stats model.matrix
#' @importFrom SummarizedExperiment assays
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
#' # (1) test for differentially abundant (DA) clusters
#' group <- factor(group_IDs, levels = c("ref", "BCRXL"))  # re-level factor to use "ref" as base level
#' res_DA <- testDA(d_clus, group)
#' topTable(res_DA, number = 6)
#' 
#' # plot top differentially abundant (DA) clusters
#' # note there is no evidence for DA in this example data set (data set is too small)
#' # plotTopDAClusters(res_DA)
#' 
#' # (2) test for differential expression of functional markers: method "diffcyt-med"
#' res_DE_med <- testDE_med(d_clus, group)
#' # topTable: use 'coef = 2' for contrast of interest (BCRXL vs. ref)
#' topTable(res_DE_med, coef = 2, number = 6)
testDE_med <- function(d_clus, group, 
                       plot = TRUE, path = ".", filename = "results_DE_diffcyt_med.pdf") {
  
  if (!is.factor(group)) group <- factor(group, levels = unique(group))
  
  # model matrix for limma
  mm <- model.matrix(~ group)
  
  # number of cells per cluster
  counts <- assays(d_clus)[["n_cells"]]
  
  # extract medians and create concatenated matrix for limma
  func_names <- names(assays(d_clus))[-match("n_cells", names(assays(d_clus)))]  # "n_cells" is last
  meds <- do.call("rbind", as.list(assays(d_clus)[func_names]))
  
  weights <- do.call("rbind", rep(list(counts), length(func_names)))
  
  # fit linear models using limma
  fit <- lmFit(meds, design = mm, weights = weights)
  
  # calculate empirical Bayes moderated tests
  fit <- eBayes(fit)
  
  if (plot) {
    pdf(file.path(path, filename), width = 7, height = 7)
    plot(fit$Amean, fit$sigma, pch = ".")
    dev.off()
  }
  
  fit
}


