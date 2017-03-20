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
#' @param d_counts \code{\link[SummarizedExperiment]{SummarizedExperiment}} object 
#'   containing cluster counts (frequencies), from \code{\link{calcCounts}}.
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
#' # generate clusters (note: using small number of clusters for demonstration purposes)
#' d_se <- generateClusters(d_se, cols_to_use = lineage_cols, xdim = 4, ydim = 4, seed = 123)
#' # plotMST(d_se)
#' 
#' # calculate cluster counts
#' d_counts <- calcCounts(d_se)
#' 
#' # calculate cluster medians
#' d_medians <- calcMedians(d_se)
#' 
#' # calculate ECDFs
#' d_ecdfs <- calcECDFs(d_se)
#' 
#' 
#' ########################################################
#' # (1) Test for differential abundance (DA) of clusters #
#' ########################################################
#' 
#' # re-level factor to use "ref" as base level
#' group <- factor(group_IDs, levels = c("ref", "BCRXL"))
#' 
#' res_DA <- testDA(d_counts, group)
#' 
#' topTable(res_DA, number = 6)
#' 
#' # plot top DA clusters
#' # (note: this is a small example data set used for demonstration purposes only; results
#' # are not biologically meaningful)
#' # plotTopDAClusters(res_DA)
#' 
testDA <- function(d_counts, group) {
  
  if (!is.factor(group)) group <- factor(group, levels = unique(group))
  
  # transform counts to square-root proportional counts per sample
  counts <- assays(d_counts)[[1]]
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


