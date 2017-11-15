#' Create contrast matrix
#' 
#' Create contrast matrix for differential testing
#' 
#' Creates a contrast matrix specifying the comparison of interest, in the correct format
#' for the differential testing functions. This can then be provided to the differential
#' testing functions, together with either a design matrix or model formula, and the data
#' object.
#' 
#' The argument \code{contrast} defines the contrast of interest for the levels of
#' \code{group_IDs}. Note that the model formula (\code{\link{createFormula}}) or design
#' matrix (\code{\link{createDesignMatrix}}) assumes that the first level of
#' \code{group_IDs} is the reference level, and the model formula or design matrix are
#' specified so that the first column of the design matrix is an intercept. In this setup,
#' the subsequent terms in the design matrix represent the difference compared to the
#' reference.
#' 
#' For example, to compare the second level against the first level for a \code{group_IDs}
#' vector of length 2, with no other columns in the design matrix, the \code{contrast}
#' argument would be \code{c(0, 1)}. To compare the second level against the first level
#' for a \code{group_IDs} vector of length 3, it would be \code{c(0, 1, 0)}. If there are
#' additional terms in the design matrix (such as block IDs), which are not of interest
#' for inference, extra zeros should be included for the corresponding columns of the
#' design matrix.
#' 
#' The contrast matrix is formatted as a matrix with a single column containing the
#' contrast of interest. To perform tests for multiple contrasts, run this function and
#' the corresponding differential testing function multiple times.
#' 
#' 
#' @param group_IDs Vector or factor of group membership labels for each sample (e.g.
#'   diseased vs. healthy, or treated vs. untreated). Vectors are converted to factors
#'   internally. The first level of the factor (e.g. healthy) or first entry of the vector
#'   will be used as the reference level for differential testing. To re-order factor
#'   levels, use \code{\link[stats]{relevel}} or \code{\link[base]{factor}}.
#' 
#' @param contrast Vector defining the contrast of interest. For example, \code{c(0, 1)}
#'   to compare the second level against the first level for a \code{group_IDs} vector of
#'   length 2, with no other columns in the design matrix. See details for additional
#'   explanation. Default is to compare the second level against the first level of
#'   \code{group_IDs}, assuming no other columns in the design matrix.
#' 
#' 
#' @return Returns a contrast matrix containing the contrast of interest, formatted as a
#'   matrix with a single column.
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
createContrast <- function(group_IDs, contrast = NULL) {
  
  if (!is.factor(group_IDs)) {
    group_IDs <- factor(group_IDs, levels = unique(group_IDs))
  }
  
  # default: compare second level against first level of 'group_IDs'
  if (is.null(contrast)) {
    contrast <- rep(0, length(levels(group_IDs)))
    contrast[2] <- 1
  }
  
  # create contrast matrix
  contrast_matrix <- matrix(contrast, ncol = 1)
  row_nms <- seq_along(contrast)
  row_nms[seq_len(nlevels(group_IDs))] <- levels(group_IDs)
  rownames(contrast_matrix) <- row_nms
  
  contrast_matrix
}


