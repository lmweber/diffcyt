#' Transform data
#' 
#' Transform data prior to clustering
#' 
#' Flow and mass cytometry data should be transformed prior to clustering. The raw data
#' follows an approximately log-normal distribution. Transforming with a log (or similar)
#' function brings the data closer to a normal distribution, which improves clustering
#' performance and allows positive and negative populations to be visualized more clearly.
#' 
#' This function implements an inverse hyperbolic sine ('arcsinh') transform with
#' adjustable 'cofactor' parameter. The arcsinh transform is widely used for CyTOF data.
#' It behaves similarly to a log transform at high values, but is approximately linear
#' near zero; so unlike the log, it can handle zeros or small negative values. The
#' cofactor parameter controls the width of the linear region. Zero values and small
#' negatives occur in CyTOF data when no ions are detected in a given channel (negatives
#' are due to background subtraction and randomization of integer count values, which are
#' performed by default by the CyTOF instrument software).
#' 
#' Recommended values for the cofactor parameter are 5 (mass cytometry, CyTOF) or 150
#' (fluorescence flow cytometry); see Bendall et al. (2011), \emph{Science}, Supplementary
#' Figure S2.
#' 
#' The transform should be applied to protein marker columns only. The
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} object created in the previous
#' step (\code{\link{prepareData}}) is assumed to contain a vector of logical entries
#' (\code{is_marker_col}) in the column meta-data, indicating which columns are marker
#' columns. (If this is not available, all columns will be transformed instead.)
#' 
#' 
#' @param d_se Input data. Assumed to be in the form of a
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}}, prepared with the function
#'   \code{\link{prepareData}}. Column meta-data is assumed to contain a vector of logical
#'   entries (\code{is_marker_col}) indicating marker columns.
#' 
#' @param cofactor Cofactor parameter for 'arcsinh' transform. Default = 5, which is
#'   appropriate for mass cytometry (CyTOF) data. For fluorescence flow cytometry, we
#'   recommend cofactor = 150 instead.
#' 
#' 
#' @return d_se Data with transform applied to protein marker columns.
#' 
#' 
#' @importFrom SummarizedExperiment assays colData 'assays<-'
#' 
#' @export
#'
#' @seealso \code{\link{testDA}}
#'
#' @examples
#' # See full examples in testing functions.
#' 
transformData <- function(d_se, cofactor = 5) {
  
  is_marker_col <- colData(d_se)$is_marker_col
  
  if (is.null(is_marker_col)) is_marker_col <- rep(TRUE, ncol(d_se[[1]]))
  
  # extract expression data
  d_ex <- assays(d_se)[[1]]
  
  # apply transform to marker columns
  d_ex[, is_marker_col] <- asinh(d_ex[, is_marker_col] / cofactor)
  
  assays(d_se)[[1]] <- d_ex
  
  d_se
}


