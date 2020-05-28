# Some checks to make sure the data is in the right format, type conversions
#
#' @importFrom magrittr %>%
#' @importFrom rlang :=
data_processing_for_imputation <- function(data, censored_variable,
                                           censoring_indicator,
                                           response = NULL, covariates = NULL,
                                           id = NULL){
  data <- dplyr::as_tibble(data)
  # transform censored variable to double
  if (is.integer(data[[censored_variable]])){
    data <- dplyr::mutate(data, 
                          !!censored_variable := as.double(data[[!!censored_variable]]))
  }
  # transform censoring indicator to integer and check that it is 0 or 1
  if (is.factor(purrr::as_vector(data[[censoring_indicator]]))){
    data <- dplyr::mutate(data,
                          !!censoring_indicator := as.integer(data[[!!censoring_indicator]])-1)
  }
  if (!is.integer(purrr::as_vector(data[[censoring_indicator]]))){
    data <- dplyr::mutate(data, 
                          !!censoring_indicator := as.integer(data[[!!censoring_indicator]]))
  }
  if (!all(data[[censoring_indicator]] %in% c(0,1))){
    rlang::abort("Not all values of 'censoring_indicator' are 0 or 1")
  }
  # transform covariate from factor to double
  if (!is.null(covariates)){
    data <- covariates_from_factor_to_numeric(data, covariates)
  }
  # sort according to censored_variable
  data <- data %>%
    dplyr::arrange(!!dplyr::sym(censored_variable)) 
  return(data)
}

# transform factor covariates to numeric
covariates_from_factor_to_numeric <- function(data, covariates){
  for (cov in covariates){
    if (is.factor(data[[cov]])){
      if (length(levels(data[[cov]])) <= 2){
        data[[cov]] <- as.numeric(data[[cov]])-1
        stopifnot(all(data[[cov]] %in% c(0,1)))
      } else{
        stop("multi-level factors not supported")
      }
    }
  }
  return(data)
}

# Extract regression coefficients beta from a regression fit, lm-like object
# 
#' @importFrom stats coef
get_betas_from_fit <- function(fit, regression_type){
  betas <- tryCatch({
    if (!identical(regression_type,"glmer")) {
      coef(fit)
    } else{
      lme4::fixef(fit)
    }
  }, error = function(cnd){lme4::fixef(fit)})
  names(betas) <- paste0("b",seq_along(betas)-1)
  return(betas)
}


# Extract regression coefficients beta from a list of regression fits
get_betas_from_multiple_fits <- function(fits_list, regression_type){
  betas_ls <- purrr::map(fits_list, ~ get_betas_from_fit(.x, regression_type))
  return(do.call(rbind,betas_ls))
}


# Extract variance of regression coefficients beta from a regression fit, 
# lm-like object
# 
#' @importFrom stats vcov
get_variance_of_betas_from_fit <- function(fit){
  var_betas <- diag(as.matrix(vcov(fit)))
  names(var_betas) <- paste0("Var(b",seq_along(var_betas)-1, ")")
  return(var_betas)
}


# Extract variances of regression coefficients beta from a list of regression fits
get_variance_of_betas_from_multiple_fits <- function(fits_list){
  var_betas_ls <- purrr::map(fits_list, ~ get_variance_of_betas_from_fit(.x))
  return(do.call(rbind,var_betas_ls))
}


# create arguments list for regression fitting
# 
# for weights to be used, must be a column in 'data' with name weights
args_for_fitting <- function(data, formula, regression_type, family = "binomial", weights = NULL){
  stopifnot(regression_type %in% c("lm","glm","glmer"))
  args = list(formula = formula, data = data)
  if (!identical(regression_type,"lm")) {
    args[["family"]] <- family
    if (!is.null(weights)) {
      args[["weights"]] <- data[[weights]]
    }
  }
  return(args)
}


# Logical. Check if the highest censored value is censored
last_is_censored <- function(data, censoring_indicator){
  return(data[dim(data)[1],censoring_indicator] == 0)
}


# set the highest censored value as observed
set_last_as_observed <- function(data, censored_variable, censoring_indicator){
  if (last_is_censored(data, censoring_indicator)){
    data[dim(data)[1],censoring_indicator] <- 1
    # repeat for duplicate values from bootstrapping
    data[data[[censored_variable]] == data[[censored_variable]][dim(data)[1]],censoring_indicator] <- 1
  }
  return(data)
}


# set the highest censored value as censored
set_last_as_censored <- function(data, censored_variable, censoring_indicator){
  if (!last_is_censored(data, censoring_indicator)){
    data[dim(data)[1],censoring_indicator] <- 0
    # repeat for duplicate values from bootstrapping
    data[data[[censored_variable]] == data[[censored_variable]][dim(data)[1]],censoring_indicator] <- 0
  }
  return(data)
}


create_two_level_factor_data <- function(data, covariates, max_levels = 4){
  stopifnot(all(covariates %in% colnames(data)))
  # loop through covariates
  for (cov in covariates) {
    if (is.factor(data[[cov]])){
      if (length(levels(data[[cov]]))<=max_levels){
        # loop through factor levels, total number is number of levels -1
        for (i in seq_along(levels(data[[cov]])[-1])){
          # new name
          cov_name <- paste0(cov, "_",levels(data[[cov]])[i+1])
          # add to data frame
          data[[cov_name]] <- factor(as.integer(data[[cov]] == levels(data[[cov]])[i+1]))
        } 
      }
    }
  }
  return(data)
}


create_two_level_factor_covariates <- function(data, formula, max_levels = 4){
  variables <- extract_variables_from_formula(formula)
  covariates <- variables$covariates
  stopifnot(all(covariates %in% colnames(data)))
  new_covariates <- c()
  # loop through covariates
  for (cov in covariates) {
    if (length(levels(data[[cov]]))<=max_levels){
      if (is.factor(data[[cov]])){
      # loop through factor levels, total number is number of levels -1
        for (i in seq_along(levels(data[[cov]])[-1])){
          # new name
          cov_name <- paste0(cov, "_",levels(data[[cov]])[i+1])
          new_covariates <- c(new_covariates, cov_name)
        }
      }
    } 
  }
  return(new_covariates)
}
