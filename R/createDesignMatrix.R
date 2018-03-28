#' Create design matrix
#' 
#' Create design matrix for model fitting
#' 
#' Creates a design matrix specifying the models to be fitted. (Alternatively,
#' \code{\link{createFormula}} can be used to generate a model formula instead of a design
#' matrix.)
#' 
#' The design matrix can then be provided to the differential testing functions, together
#' with the data object and contrast matrix.
#' 
#' The \code{sample_info} input (which was also previously provided to
#' \code{\link{prepareData}}) should be a data frame containing all factors and covariates
#' of interest. For example, depending on the experimental design, this may include the
#' following columns:
#' 
#' \itemize{
#' \item group IDs (e.g. groups for differential testing)
#' \item block IDs (e.g. patient IDs in a paired design)
#' \item batch IDs (batch effects)
#' \item continuous covariates
#' }
#' 
#' The logical vector \code{cols_include} specifies which columns in \code{sample_info} to
#' include in the design matrix. (For example, there may be an additional column of sample
#' IDs, which should not be included.)
#' 
#' Columns of indicator variables (e.g. group IDs, block IDs, and batch IDs) must be
#' formatted as factors (otherwise they will be treated as numeric values). The indicator
#' columns will be expanded into the design matrix format. The names for each parameter
#' are taken from the column names of \code{sample_info}.
#' 
#' All factors provided here will be included as fixed effect terms in the design matrix.
#' Alternatively, to use random effects for some factors (e.g. for block IDs), see
#' \code{\link{createFormula}}; or, depending on the method used, provide them directly to
#' the differential testing function (\code{\link{testDA_limma}} and
#' \code{\link{testDS_limma}}).
#' 
#' 
#' @param sample_info Data frame of sample information (which was also previously provided
#'   to \code{\link{prepareData}}). This should be a data frame containing all factors and
#'   covariates of interest; e.g. group IDs, block IDs, batch IDs, and continuous
#'   covariates.
#' 
#' @param cols_include (Logical) Columns of \code{sample_info} to include in the design
#'   matrix. Default = all columns.
#' 
#' 
#' @return \code{design}: Returns a design matrix (numeric matrix), with one row per
#'   sample, and one column per model parameter.
#' 
#' 
#' @importFrom stats as.formula model.matrix
#' @importFrom methods is
#' 
#' @export
#' 
#' @examples
#' # For a full workflow example demonstrating the use of each function in the 'diffcyt'
#' # pipeline, see the package vignette.
#' 
#' # Example: simple design matrix
#' sample_info <- data.frame(
#'   sample_IDs = paste0("sample", 1:4), 
#'   group_IDs = factor(c("group1", "group1", "group2", "group2"))
#' )
#' createDesignMatrix(sample_info, cols_include = 2)
#' 
#' # Example: more complex design matrix: patient IDs and batch IDs
#' sample_info <- data.frame(
#'   sample_IDs = paste0("sample", 1:8), 
#'   group_IDs = factor(rep(paste0("group", 1:2), each = 4)), 
#'   patient_IDs = factor(rep(paste0("patient", 1:4), 2)), 
#'   batch_IDs = factor(rep(paste0("batch", 1:2), 4))
#' )
#' createDesignMatrix(sample_info, cols_include = 2:4)
#' 
#' # Example: more complex design matrix: continuous covariate
#' sample_info <- data.frame(
#'   sample_IDs = paste0("sample", 1:4), 
#'   group_IDs = factor(c("group1", "group1", "group2", "group2")), 
#'   age = c(52, 35, 71, 60)
#' )
#' createDesignMatrix(sample_info, cols_include = 2:3)
#' 
createDesignMatrix <- function(sample_info, cols_include = NULL) {
  
  stopifnot(is.data.frame(sample_info))
  
  if (is.null(cols_include)) cols_include <- seq_len(nrow(sample_info))
  
  # create design matrix
  terms <- colnames(sample_info)[cols_include]
  
  formula <- as.formula(paste("~", paste(terms, collapse = " + ")))
  
  design <- model.matrix(formula, data = sample_info)
  
  design
}


