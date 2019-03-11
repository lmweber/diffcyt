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
#' zero. In many cases, this will simply be a vector of zeros and a single entry equal to
#' one; this will test whether a single parameter is equal to zero (e.g. c(0, 1, 0, 0,
#' 0)).
#' 
#' If a design matrix has been used, the entries of \code{contrast} correspond to the
#' columns of the design matrix; and the length of \code{contrast} equals the number of
#' columns in the design matrix. If a model formula has been used, the entries correspond
#' to the levels of the fixed effect terms; and the length equals the number of levels of
#' the fixed effect terms.
#' 
#' The contrast matrix is formatted as a matrix with a single column containing the
#' contrast of interest. To perform tests for multiple contrasts, run this function and
#' the corresponding differential testing function multiple times.
#' 
#' 
#' @param contrast Vector defining the contrast of interest. This should be a numeric
#'   vector specifying the combination of model parameters to test whether they are equal
#'   to zero. The entries correspond to the columns of the design matrix, or the levels of
#'   the fixed effect terms in the model formula. For example, using a design matrix:
#'   \code{c(0, 1, 0, 0, 0)} to test whether a single parameter corresponding to the
#'   second column in the design matrix is equal to zero.
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
#' # For a complete workflow example demonstrating each step in the 'diffcyt' pipeline, 
#' # see the package vignette.
#' 
#' # Example: contrast matrix
#' createContrast(c(0, 1, 0, 0, 0))
#' 
createContrast <- function(contrast) {
  
  contrast_matrix <- matrix(contrast, ncol = 1)
  
  contrast_matrix
}


