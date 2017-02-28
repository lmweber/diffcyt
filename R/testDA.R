#' Test clusters for differential abundance
#' 
#' Calculate empirical Bayes moderated tests for differential abundance of clusters
#' 
#' Calculates tests for differential abundance (differential frequencies) of clusters, 
#' using empirical Bayes moderation of cluster-wise variances to improve power.
#' 
#' We use the \code{\link[limma]{limma}} package (Ritchie et al. 2015, \emph{Nucleic
#' Acids Research}) to calculate the empirical Bayes moderated tests. Empirical Bayes
#' methods improve statistical power by sharing information on variability (i.e. variance
#' across samples for a single cluster) across clusters.
#' 
#' This function should be run after generating clusters with
#' \code{\link{generateClusters}}.
#' 
#' 
#' @param d_transf Transformed input data from the previous steps. This should be in the
#'   form of a \code{\link[flowCore]{flowSet}} object from the
#'   \code{\link[flowCore]{flowCore}} package.
#' 
#' @param group Factor containing group membership for each sample (for example, diseased
#'   vs. healthy). This is required to fit linear models and calculate statistical tests 
#'   for each cluster. [to do: update this to extract it automatically from the flowSet 
#'   object]
#' 
#' @param tbl_prop Table of cluster proportions (rows = clusters, columns = samples).
#'   Calculated by function \code{\link{calculateFreq}}. [to do: update to save tbl_freq
#'   and tbl_prop directly in an object, possibly as an additional sheet in a Summarized
#'   Experiment]
#' 
#' 
#' @return Returns \code{\link[limma]{limma}} fitted model objects for each cluster,
#'   including p-values from empirical Bayes moderated tests. The p-values can be
#'   accessed with \code{\link[limma]{topTable}} from the \code{limma} package.
#' 
#' 
#' @importFrom limma lmFit eBayes
#' @importFrom stats model.matrix
#' @importFrom methods as is
#' 
#' @export
#' 
#' @examples
#' # need to create a small example data set for examples
testDA <- function(d_transf, group, tbl_prop) {
  
  # table of proportions
  tbl_prop <- t(t(tbl_freq) / colSums(tbl_freq))  # transpose because vector wraps by column
  stopifnot(all(colSums(tbl_prop) == 1))
  
  stopifnot(is(group, "factor"))
  
  sample_IDs <- rownames(flowCore::phenoData(d_transf))
  
  # model matrix for limma
  mm <- model.matrix(~ group)
  
  # fit linear models on transformed proportions (square root)
  f1 <- limma::lmFit(sqrt(tbl_prop), design = mm)
  
  # add proportions (expressed as percentages) to fitted model object
  p <- matrix(100 * as.numeric(tbl_prop), ncol = length(sample_IDs))  # needs to be simple matrix
  colnames(p) <- colnames(tbl_prop)
  rownames(p) <- rownames(tbl_prop)
  f1$genes <- p
  
  # calculate empirical Bayes moderated tests
  f1 <- limma::eBayes(f1)
  
  f1
}


