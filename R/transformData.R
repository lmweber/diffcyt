#' Transform data
#' 
#' Transform data prior to clustering
#' 
#' Flow and mass cytometry data should be transformed prior to clustering. This function 
#' implements an inverse hyperbolic sine ('arcsinh') transform with adjustable 'cofactor'
#' parameter. Recommended values for the cofactor are 5 (mass cytometry, CyTOF) or 150
#' (fluorescence flow cytometry); see Bendall et al. (2011), \emph{Science}, Supplementary
#' Figure S2.
#' 
#' The transform will be applied to protein marker columns only. The 
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} object created in the previous
#' step (\code{\link{prepareData}}) is assumed to contain a vector of logical entries 
#' (\code{is_marker}) in the column meta-data, indicating which columns are marker 
#' columns; otherwise all columns are transformed.
#' 
#' 
#' @param d_se Input data. Assumed to be in the form of a 
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}}, prepared with the function 
#'   \code{\link{prepareData}}. Column meta-data is assumed to contain a vector of logical
#'   entries (\code{is_marker}) indicating marker columns.
#' 
#' @param cofactor Cofactor for 'arcsinh' transform. Default = 5, which is appropriate for
#'   mass cytometry (CyTOF) data. For fluorescence flow cytometry, we recommend cofactor =
#'   150 instead.
#' 
#' 
#' @return d_se Data with transform applied to protein marker columns.
#' 
#' 
#' @importFrom SummarizedExperiment assays colData 'assays<-'
#' 
#' @export
#'
#' @seealso \code{\link{testDA}} \code{\link{testDE_med}} \code{\link{testDE_FDA}}
#'
#' @examples
#' # See full examples in testing functions: testDA, testDE_med, testDE_FDA
#' 
transformData <- function(d_se, cofactor = 5) {
  
  is_marker <- colData(d_se)$is_marker
  
  if (is.null(marker_cols)) marker_cols <- 1:ncol(d_se[[1]])
  
  # extract expression data
  d_ex <- assays(d_se)[[1]]
  
  # transform marker columns
  d_ex[, is_marker] <- asinh(d_ex[, is_marker] / cofactor)
  
  assays(d_se)[[1]] <- d_ex
  
  d_se
}


