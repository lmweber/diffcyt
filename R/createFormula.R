#' Create model formula and corresponding data frame of variables
#' 
#' Create model formula and corresponding data frame of variables for model fitting
#' 
#' Creates a model formula and corresponding data frame of variables specifying the models
#' to be fitted. (Alternatively, \code{\link{createDesignMatrix}} can be used to generate
#' a design matrix instead of a model formula.)
#' 
#' The output is a list containing the model formula and corresponding data frame of
#' variables (one column per formula term). These can then be provided to differential
#' testing functions that require a model formula, together with the main data object and
#' contrast matrix.
#' 
#' The \code{experiment_info} input (which was also previously provided to
#' \code{\link{prepareData}}) should be a data frame containing all factors and covariates
#' of interest. For example, depending on the experimental design, this may include the
#' following columns:
#' 
#' \itemize{
#' \item group IDs (e.g. groups for differential testing)
#' \item block IDs (e.g. patient IDs in a paired design; these may be included as either
#' fixed effect or random effects)
#' \item batch IDs (batch effects)
#' \item continuous covariates
#' \item sample IDs (e.g. to include random intercept terms for each sample, to account
#' for overdispersion typically seen in high-dimensional cytometry data; this is known as
#' an 'observation-level random effect' (OLRE); see see Nowicka et al., 2017,
#' \emph{F1000Research} for more details)
#' }
#' 
#' The arguments \code{cols_fixed} and \code{cols_random} specify the columns in
#' \code{experiment_info} to include as fixed effect terms and random intercept terms
#' respectively. These can be provided as character vectors of column names, numeric
#' vectors of column indices, or logical vectors. The names for each formula term are
#' taken from the column names of \code{experiment_info}.
#' 
#' Note that for some methods, random effect terms (e.g. for block IDs) must be provided
#' directly to the differential testing function instead (\code{\link{testDA_voom}} and
#' \code{\link{testDS_limma}}).
#' 
#' If there are no random effect terms, it will usually be simpler to use a design matrix
#' instead of a model formula; see \code{\link{createDesignMatrix}}.
#' 
#' 
#' @param experiment_info \code{data.frame}, \code{DataFrame}, or \code{tbl_df} of
#'   experiment information (which was also previously provided to
#'   \code{\link{prepareData}}). This should be a data frame containing all factors and
#'   covariates of interest; e.g. group IDs, block IDs, batch IDs, and continuous
#'   covariates.
#' 
#' @param cols_fixed Argument specifying columns of \code{experiment_info} to include as
#'   fixed effect terms in the model formula. This can be provided as a character vector
#'   of column names, a numeric vector of column indices, or a logical vector.
#' 
#' @param cols_random Argument specifying columns of \code{experiment_info} to include as
#'   random intercept terms in the model formula. This can be provided as a character
#'   vector of column names, a numeric vector of column indices, or a logical vector.
#'   Default = none.
#' 
#' 
#' @return \code{formula}: Returns a list with three elements:
#' \itemize{
#' \item \code{formula}: model formula
#' \item \code{data}: data frame of variables corresponding to the model formula
#' \item \code{random_terms}: TRUE if model formula contains any random effect terms
#' }
#' 
#' 
#' @importFrom stats as.formula
#' @importFrom methods is
#' 
#' @export
#' 
#' @examples
#' # For a complete workflow example demonstrating each step in the 'diffcyt' pipeline, 
#' # see the package vignette.
#' 
#' # Example: model formula
#' experiment_info <- data.frame(
#'   sample_id = factor(paste0("sample", 1:8)), 
#'   group_id = factor(rep(paste0("group", 1:2), each = 4)), 
#'   patient_id = factor(rep(paste0("patient", 1:4), 2)), 
#'   stringsAsFactors = FALSE
#' )
#' createFormula(experiment_info, cols_fixed = "group_id", cols_random = c("sample_id", "patient_id"))
#' 
createFormula <- function(experiment_info, cols_fixed = NULL, cols_random = NULL) {
  
  stopifnot(any(class(experiment_info) %in% c("data.frame", "DataFrame", "tbl_df", "tbl")))
  experiment_info <- as.data.frame(experiment_info)
  
  # create formula
  
  # fixed effects
  if (is.character(cols_fixed)) {
    stopifnot(all(cols_fixed %in% colnames(experiment_info)))
    terms_fixed <- cols_fixed
  } else if (is.numeric(cols_fixed) | is.logical(cols_fixed)) {
    terms_fixed <- colnames(experiment_info)[cols_fixed]
  }
  RHS_fixed <- paste(terms_fixed, collapse = " + ")
  
  # random effects
  if (!is.null(cols_random)) {
    if (is.character(cols_random)) {
      stopifnot(all(cols_random %in% colnames(experiment_info)))
      terms_random <- cols_random
    } else if (is.numeric(cols_random) | is.logical(cols_random)) {
      terms_random <- colnames(experiment_info)[cols_random]
    }
    RHS_random <- paste(paste0("(1 | ", terms_random, ")"), collapse = " + ")
    random_terms <- TRUE
  } else if (is.null(cols_random)) {
    terms_random <- NULL
    RHS_random <- NULL
    random_terms <- FALSE
  }
  
  # combine
  RHS_combined <- paste(c(RHS_fixed, RHS_random), collapse = " + ")
  
  formula <- as.formula(paste(c("y", RHS_combined), collapse = " ~ "))
  
  # create data frame of corresponding variables (in correct order)
  cols <- c(terms_fixed, terms_random)
  data <- experiment_info[, cols, drop = FALSE]
  
  # return as list
  list(formula = formula, data = data, random_terms = random_terms)
}


