#' check if 'Surv()' formula is valid
#'
#' @param formula formula containing indication 'Surv()'.
#'  See \code{\link{testDA_censoredGLMM}}
#' @param throw_error should an error be thrown if 'formula' is incorrect.
#'  default = 'TRUE'.
#' @return logical
#' @export
is_valid_censored_formula <- function(formula, throw_error = TRUE){
  tmpfor_char <- as.character(formula)
  if (length(tmpfor_char) != 3) {
    if (throw_error) stop("no response specified in formula",call. = FALSE)
    return(FALSE)
  }
  which_surv <- stringr::str_detect(tmpfor_char,"Surv")
  if (which_surv[2]) {
    if (throw_error) stop("censored response is not allowed in formula",call. = FALSE)
    return(FALSE)
  }
  if (!any(which_surv)) {
    if (throw_error) stop("no censored covariate specified in formula",call. = FALSE)
    return(FALSE)
  }
  return(TRUE)
}


#' extract individuall components of 'Surv()' formula
#'
#' @param formula formula containing indication 'Surv()'.
#'  See \code{\link{testDA_censoredGLMM}}
#' @return list with elements 'censored_variable', 'censoring_indicator',
#'  'response', 'covariates', 'random_covariates'. each a character vector.
#' @export
extract_variables_from_formula <- function(formula){
  is_valid_censored_formula(formula = formula, throw_error = TRUE)
  tmpfor_char <- as.character(formula)
  which_surv <- stringr::str_detect(tmpfor_char,"Surv")
  right_side <- stringr::str_split(tmpfor_char[which_surv],"[:space:]*\\+[:space:]*")[[1]]
  which_surv <- stringr::str_detect(right_side,"Surv")
  tmpfor_char_surv <- stringr::str_extract(right_side[which_surv],"\\([[:graph:][:space:]]+(?=\\))")
  tmpfor_char_surv <- stringr::str_replace(tmpfor_char_surv,"\\(","")
  tmpfor_char_surv <- stringr::str_split(tmpfor_char_surv,"[:space:]*,[:space:]*")[[1]]
  censored_variable <- tmpfor_char_surv[1]
  censoring_indicator <- tmpfor_char_surv[2]
  covariates <- stringr::str_extract(right_side[!which_surv],"([:alpha:]+[:alnum:]*[_.]*)+([:alpha:]*[:alnum:]*[_.]*)*")
  covariates <- covariates[covariates != "" & !is.na(covariates)]
  is_random_term <- stringr::str_detect(right_side[!which_surv], "\\|")
  random_covariates <- covariates[is_random_term]
  if (identical(random_covariates,character(0))) {
    random_covariates <- NULL
  }
  covariates <- covariates[!is_random_term]
  if (identical(covariates,character(0))) {
    covariates <- NULL
  }
  response <- tmpfor_char[2]
  return(list(censored_variable = censored_variable, censoring_indicator = censoring_indicator,
              response = response, covariates = covariates, random_covariates = random_covariates))
}

#' create formula from formula containing 'Surv()'
#'
#' @param formula formula containing indication 'Surv()'.
#'  See \code{\link{testDA_censoredGLMM}}
#' @return formula containing no 'Surv()' terms.
#' @export
create_glmm_formula <- function(formula){
  extracted_vars <- extract_variables_from_formula(formula)
  outformula <- stringr::str_replace(formula, "Surv\\([[:alnum:] ,.]*\\)",extracted_vars$censored_variable)
  as.formula(paste0(outformula[2],outformula[1],outformula[3]))
}


