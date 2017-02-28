#' Test clusters for differential abundance
#' 
#' Calculate empirical Bayes moderated tests for differential abundance of clusters
#' 
#' Calculates tests for differential abundance (differential frequencies) of clusters, 
#' using empirical Bayes moderation of cluster-wise variances to improve power.
#' 
#' We use the \code{\link[limma]{limma}} package (Ritchie et al. 2015, \emph{Nucleic Acids
#' Research}) to calculate the empirical Bayes moderated tests. Empirical Bayes methods 
#' improve statistical power by sharing information on variability (i.e. variance across 
#' samples for a single cluster) across clusters.
#' 
#' Note that an additional transformation is included: the tests for differential
#' abundance are calculated using the square root of proportional counts per sample,
#' instead of raw counts.
#' 
#' 
#' @param d_clus \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'   containing cluster medians (median expression of functional markers) and cluster
#'   frequencies (number of cells), from \code{\link{calcMediansAndFreq}}.
#' 
#' @param group Factor containing group membership for each sample (for example, diseased
#'   vs. healthy), for differential comparisons and statistical tests.
#' 
#' 
#' @return Returns fitted model objects from \code{\link[limma]{limma}} for each cluster, 
#'   including p-values from the empirical Bayes moderated tests for differential
#'   abundance. The p-values can be accessed with the function
#'   \code{\link[limma]{topTable}} from the \code{limma} package.
#' 
#' 
#' @importFrom SummarizedExperiment assays
#' @importFrom limma lmFit eBayes
#' @importFrom stats model.matrix
#' @importFrom methods as is
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
testDA <- function(d_clus, group) {
  
  if (!is.factor(group)) group <- factor(group, levels = unique(group))
  
  # transform counts to square-root proportional counts per sample
  counts <- assays(d_clus)[["n_cells"]]
  prop <- t(t(counts) / colSums(counts))  # transpose because colSums vector wraps by column
  stopifnot(all(colSums(prop) == 1))
  sqrt_prop <- sqrt(prop)
  
  # model matrix for limma
  mm <- model.matrix(~ group)
  
  # fit linear models using limma
  fit <- lmFit(sqrt_prop, design = mm)
  
  # include proportions (expressed as percentages) in fitted model object
  p <- matrix(100 * as.numeric(prop), ncol = length(group))  # needs to be simple matrix
  colnames(p) <- colnames(prop)
  rownames(p) <- rownames(prop)
  fit$genes <- p
  
  # calculate empirical Bayes moderated tests
  fit <- eBayes(fit)
  
  fit
}


