#' Test for differential expression of functional markers (method: "diffcyt-med")
#' 
#' Calculate tests for differential expression of functional markers within clusters 
#' (method: "diffcyt-med"), using empirical Bayes moderation of variances to improve
#' power.
#' 
#' The "diffcyt-med" methodology uses median functional marker expression to characterize
#' the signal within each cluster. The differential expression tests compare the median
#' functional marker expression within clusters between samples in the two groups (e.g.
#' diseased vs. healthy).
#' 
#' The number of cells per cluster per sample is used to weight the tests, representing
#' the uncertainty in calculating each median value.
#' 
#' Filtering: Clusters are kept for testing if there are at least \code{min_samples}
#' samples in each condition with at least \code{min_cells} cells. Filtered clusters are
#' removed from testing for all functional markers.
#' 
#' Empirical Bayes methods are used to share information on variability across clusters, 
#' improving power. We use the \code{\link[limma]{limma}} package (Ritchie et al. 2015, 
#' \emph{Nucleic Acids Research}) to calculate the empirical Bayes moderated tests.
#' 
#' Alternative methodologies for testing for differential expression of functional markers
#' are available in the functions \code{testDE-FDA} and \code{testDE-KS}.
#' 
#' 
#' @param d_medians \code{\link[SummarizedExperiment]{SummarizedExperiment}} object 
#'   containing cluster medians (median expression of functional markers), from
#'   \code{\link{calcMedians}}.
#' 
#' @param d_counts \code{\link[SummarizedExperiment]{SummarizedExperiment}} object 
#'   containing cluster counts (frequencies), from \code{\link{calcCounts}}.
#' 
#' @param group Factor containing group membership for each sample (for example, diseased
#'   vs. healthy), for differential comparisons and statistical tests.
#' 
#' @param min_cells Filtering parameter. Default = 6. Clusters are kept if there are at
#'   least \code{min_samples} samples in each condition with at least \code{min_cells}
#'   cells. Filtered clusters are removed from testing for all functional markers.
#' 
#' @param min_samples Filtering parameter. Default = 2. Clusters are kept if there are at 
#'   least \code{min_samples} samples in each condition with at least \code{min_cells} 
#'   cells. Filtered clusters are removed from testing for all functional markers.
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
#' #########################################################################################
#' # (2) Test for differential expression (DE) of functional markers: method "diffcyt-med" #
#' #########################################################################################
#' 
#' # re-level factor to use "ref" as base level
#' group <- factor(group_IDs, levels = c("ref", "BCRXL"))
#' 
#' res_DE_med <- testDE_med(d_medians, d_counts, group)
#' 
#' # topTable: use 'coef = 2' for contrast of interest (BCRXL vs. ref)
#' # (note: this is a small example data set used for demonstration purposes only; results
#' # are not biologically meaningful)
#' topTable(res_DE_med, coef = 2, number = 6)
#' 
testDE_med <- function(d_medians, d_counts, group, min_cells = 6, min_samples = 2, 
                       plot = TRUE, path = ".", filename = "results_DE_diffcyt_med.pdf") {
  
  if (!is.factor(group)) group <- factor(group, levels = unique(group))
  
  # model matrix for limma
  mm <- model.matrix(~ group)
  
  # number of cells per cluster
  counts <- assays(d_counts)[[1]]
  
  # filtering
  grp <- group == levels(group)[1]
  tf <- counts >= min_cells
  ix_keep <- (rowSums(tf[, grp]) > min_samples) & (rowSums(tf[, !grp]) > min_samples)
  
  counts <- counts[ix_keep, ]
  
  # extract medians and create concatenated matrix for limma
  func_names <- names(assays(d_medians))
  meds <- do.call("rbind", {
    lapply(as.list(assays(d_medians)[func_names]), function(a) a[ix_keep, ])
  })
  # rownames: functional marker names and cluster labels
  rownames(meds) <- paste(rep(func_names, each = nrow(counts)), rownames(meds), sep = "_")
  
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
  
  # [to do: include IHW and topTable; give same output as testDE_FDA()]
  
  fit
}


