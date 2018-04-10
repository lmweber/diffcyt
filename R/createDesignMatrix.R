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
#' The \code{experiment_info} input (which was also previously provided to
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
#' The logical vector \code{cols_design} specifies which columns in \code{experiment_info}
#' to include in the design matrix. (For example, there may be an additional column of
#' sample IDs, which should not be included.)
#' 
#' Columns of indicator variables (e.g. group IDs, block IDs, and batch IDs) must be
#' formatted as factors (otherwise they will be treated as numeric values). The indicator
#' columns will be expanded into the design matrix format. The names for each parameter
#' are taken from the column names of \code{experiment_info}.
#' 
#' All factors provided here will be included as fixed effect terms in the design matrix.
#' Alternatively, to use random effects for some factors (e.g. for block IDs), see
#' \code{\link{createFormula}}; or, depending on the method used, provide them directly to
#' the differential testing function (\code{\link{testDA_voom}} and
#' \code{\link{testDS_limma}}).
#' 
#' 
#' @param experiment_info \code{data.frame} or \code{DataFrame} of experiment information
#'   (which was also previously provided to \code{\link{prepareData}}). This should be a
#'   data frame containing all factors and covariates of interest; e.g. group IDs, block
#'   IDs, batch IDs, and continuous covariates.
#' 
#' @param cols_design (Logical) Columns of \code{experiment_info} to include in the design
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
#' experiment_info <- data.frame(
#'   sample_id = factor(paste0("sample", 1:4)), 
#'   group_id = factor(c("group1", "group1", "group2", "group2")), 
#'   stringsAsFactors = FALSE
#' )
#' createDesignMatrix(experiment_info, cols_design = 2)
#' 
#' # Example: more complex design matrix: patient IDs and batch IDs
#' experiment_info <- data.frame(
#'   sample_id = factor(paste0("sample", 1:8)), 
#'   group_id = factor(rep(paste0("group", 1:2), each = 4)), 
#'   patient_id = factor(rep(paste0("patient", 1:4), 2)), 
#'   batch_id = factor(rep(paste0("batch", 1:2), 4)), 
#'   stringsAsFactors = FALSE
#' )
#' createDesignMatrix(experiment_info, cols_design = 2:4)
#' 
#' # Example: more complex design matrix: continuous covariate
#' experiment_info <- data.frame(
#'   sample_id = factor(paste0("sample", 1:4)), 
#'   group_id = factor(c("group1", "group1", "group2", "group2")), 
#'   age = c(52, 35, 71, 60), 
#'   stringsAsFactors = FALSE
#' )
#' createDesignMatrix(experiment_info, cols_design = 2:3)
#' 
createDesignMatrix <- function(experiment_info, cols_design = NULL) {
  
  stopifnot(class(experiment_info) %in% c("data.frame", "DataFrame"))
  if (class(experiment_info) == "DataFrame") {
    experiment_info <- as.data.frame(experiment_info)
  }
  
  if (is.null(cols_design)) cols_design <- seq_len(nrow(experiment_info))
  
  # create design matrix
  terms <- colnames(experiment_info)[cols_design]
  
  formula <- as.formula(paste("~", paste(terms, collapse = " + ")))
  
  design <- model.matrix(formula, data = experiment_info)
  
  design
}


