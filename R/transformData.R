#' Transform data
#' 
#' Transform data prior to clustering
#' 
#' Flow and mass cytometry data should be transformed prior to clustering. This function 
#' implements an 'arcsinh' transform with adjustable 'cofactor' parameter. Recommended
#' values for the cofactor are 5 (mass cytometry, CyTOF) or 150 (fluorescence flow
#' cytometry); see Bendall et al. (2011), \emph{Science}, Supplementary Figure S2.
#' 
#' 
#' @param d_input Input data. Assumed to be in the form of a 
#'   \code{\link[flowCore]{flowSet}} object. The input data format can be checked using 
#'   the function \code{\link{checkInputData}}.
#' 
#' @param cofactor Cofactor for 'arcsinh' transform. Default = 5, which is appropriate 
#'   for mass cytometry (CyTOF) data. For fluorescence flow cytometry, we recommend 
#'   cofactor = 150 instead.
#' 
#' @param marker_cols Indices of protein marker columns. The transform will be applied to
#'   these columns only. Default = NULL, which transforms all columns.
#' 
#' 
#' @return d_transf Data with transform applied to protein marker columns.
#' 
#' 
#' @importFrom flowCore exprs
#' @importFrom methods as
#' 
#' @export
#'
#' @examples
#' # need to create a small example data set for examples
transformData <- function(d_input, 
                          cofactor = 5, marker_cols = NULL) {
  
  # convert flowSet to list (easier to apply transform)
  d_input <- as(d_input, "list")
  
  if (is.null(marker_cols)) marker_cols <- 1:ncol(d_input[[1]])
  
  d_transf <- lapply(d_input, function(d) {
    e <- flowCore::exprs(d)
    e[, marker_cols] <- asinh(e[, marker_cols] / cofactor)
    flowCore::exprs(d) <- e
    d
  })
  
  # convert back to flowSet
  d_transf <- as(d_transf, "flowSet")
  checkInputData(d_transf)
  
  d_transf
}


