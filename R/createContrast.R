#' Create contrast matrix
#' 
#' Create contrast matrix for differential testing
#' 
#' Creates a contrast matrix specifying the comparison of interest, in the correct format
#' for the differential testing functions. This can then be provided to the differential
#' testing functions, together with either a design matrix or model formula, and the data
#' object.
#' 
#' The argument \code{contrast} defines the contrast of interest. This should be a numeric
#' vector specifying the combination of model parameters to test whether they are equal to
#' zero. In most cases, this will simply be a vectors of zeros and a single entry equal to
#' one; this specifies a single parameter to test whether it is equal to zero (e.g. c(0,
#' 1, 0, 0, 0)).
#' 
#' If a design matrix has been used, the length of \code{contrast} should be equal to the
#' number of columns in the design matrix. If a model formula has been used, the length of
#' \code{contrast} should equal the number of fixed effect terms.
#' 
#' The contrast matrix is formatted as a matrix with a single column containing the
#' contrast of interest. To perform tests for multiple contrasts, run this function and
#' the corresponding differential testing function multiple times.
#' 
#' taken into account during inference
#' This allows their effects to be estimated during model fitting, and taken into account
#' during inference on the group ID parameters of interest.
#' 
#' 
#' @param contrast Vector defining the contrast of interest. This should be a numeric
#'   vector specifying the combination of model parameters to test whether they are equal
#'   to zero. For example, \code{c(0, 1, 0, 0, 0)} to test whether a single parameter (the
#'   second column in the design matrix) is equal to zero.
#' 
#' 
#' @return \code{contrast}: Returns a contrast matrix containing the contrast of interest,
#'   formatted as a matrix with a single column.
#' 
#' 
#' @importFrom methods is
#' 
#' @export
#' 
#' @examples
#' # A full workflow example demonstrating the use of each function in the 'diffcyt'
#' # pipeline on an experimental data set is available in the package vignette.
#' 
createContrast <- function(contrast) {
  
  contrast_matrix <- matrix(contrast, ncol = 1)
  
  contrast_matrix
}


