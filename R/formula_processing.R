# check if 'Surv()' formula is valid
#
# @param formula formula containing indication 'Surv()'.
#  See \code{\link{testDA_censoredGLMM}}
# @param throw_error should an error be thrown if 'formula' is incorrect.
#  default = 'TRUE'.
# @return logical
# @export
# @examples 
# # correct formula:
# cens_formula <- formula(Y ~ Surv(X,I) + Z)
# is_valid_censored_formula(cens_formula)
# 
# # wrong formulas:
# wrong_cens_formula_1 <- formula( ~ Surv(X,I) + Z) # no response
# is_valid_censored_formula(wrong_cens_formula_1)
# wrong_cens_formula_2 <- formula(Y ~ X + Z) # no censored variable
# is_valid_censored_formula(wrong_cens_formula_2)
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
  how_many_surv <- stringr::str_count(tmpfor_char[3],"Surv")
  if (how_many_surv != 1){
    if (throw_error) stop("currently only one censored covariate is supported",call. = FALSE)
    return(FALSE)
  }
  tmp_str_1 <- stringr::str_locate(tmpfor_char[3], "Surv\\([[:alnum:]_, ]+\\)")
  tmp_str_2 <- stringr::str_sub(tmpfor_char[3], tmp_str_1[1]+5, tmp_str_1[2]-1)
  tmp_str_3 <- stringr::str_split(tmp_str_2, "[ ]*,[ ]*")[[1]]
  tmp_str_4 <- tmp_str_3[tmp_str_3 != ""]
  if (length(tmp_str_4) < 2){
    if (throw_error) stop("censored variable, censoring indicator or ',' is 
                          missing inside 'Surv(.,.)",call. = FALSE)
    return(FALSE)
  } else if (length(tmp_str_4) > 2){
    if (throw_error) stop("invalid number or variables inside 'Surv(.,.)",call. = FALSE)
    return(FALSE)
  }
  return(TRUE)
}


# extract individuall components of 'Surv()' formula
#
# @param formula formula containing indication 'Surv()'.
#  See \code{\link{testDA_censoredGLMM}}
# @return list with elements 'censored_variable', 'censoring_indicator',
#  'response', 'covariates', 'random_covariates'. each a character vector.
# @export
# @examples 
#   cens_formula <- formula(Y ~ Surv(X,I) + Z)
#   extract_variables_from_formula(cens_formula)
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

# create formula from formula containing 'Surv()'
#
# @param formula formula containing indication 'Surv()'.
#  See \code{\link{testDA_censoredGLMM}}
# @return formula containing no 'Surv()' terms.
# @export
# @examples 
#   cens_formula <- formula(Y ~ Surv(X,I) + Z)
#   create_glmm_formula(cens_formula)
create_glmm_formula <- function(formula){
  extracted_vars <- extract_variables_from_formula(formula)
  outformula <- stringr::str_replace(formula, "Surv\\([[:alnum:] ,._]*\\)",extracted_vars$censored_variable)
  as.formula(paste0(outformula[2],outformula[1],outformula[3]))
}


