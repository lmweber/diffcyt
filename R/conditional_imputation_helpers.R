#' @importFrom magrittr %>%
#' @importFrom dplyr .data
data_processing_for_imputation <- function(data, censored_variable,
                                           censoring_indicator,
                                           response = NULL, covariates = NULL,
                                           id = NULL){
  # transform censoring indicator to integer and check that it is 0 or 1
  data <- dplyr::as_tibble(data)
  if (is.integer(data[[censored_variable]])){
    data <- dplyr::mutate(data, 
                          !!censored_variable := as.double(data[[!!censored_variable]]))
  }
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
  if (!is.null(covariates)){
    data <- covariates_from_factor_to_numeric(data, covariates)
  }
  # sort according to censored_variable
  data <- data %>%
    dplyr::arrange(!!dplyr::sym(censored_variable)) %>%
    # rank, for duplicates
    dplyr::mutate(rank = dplyr::dense_rank(!!dplyr::sym(censored_variable)))
  return(data)
}

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

get_betas_from_multiple_fits <- function(fits_list, regression_type){
  betas_ls <- purrr::map(fits_list, ~ get_betas_from_fit(.x, regression_type))
  return(do.call(rbind,betas_ls))
}
#' @importFrom stats vcov
get_variance_of_betas_from_fit <- function(fit){
  var_betas <- diag(as.matrix(vcov(fit)))
  names(var_betas) <- paste0("Var(b",seq_along(var_betas)-1, ")")
  return(var_betas)
}
get_variance_of_betas_from_multiple_fits <- function(fits_list){
  var_betas_ls <- purrr::map(fits_list, ~ get_variance_of_betas_from_fit(.x))
  return(do.call(rbind,var_betas_ls))
}

# for weights to be used, must be a column in 'data' with name weights
args_for_fitting <- function(data, formula, regression_type, family = "binomial"){
  stopifnot(regression_type %in% c("lm","glm","glmer"))
  args = list(formula = formula, data = data)
  if (!identical(regression_type,"lm")) {
    args[["family"]] <- family
    if (!is.null(data$weights)) {
      args[["weights"]] <- data$weights
    }
  }
  return(args)
}


last_is_censored <- function(data, censoring_indicator){
  return(data[dim(data)[1],censoring_indicator] == 0)
}

set_last_as_observed <- function(data, censored_variable, censoring_indicator){
  if (last_is_censored(data, censoring_indicator)){
    data[dim(data)[1],censoring_indicator] <- 1
    # repeat for duplicate values from bootstrapping
    data[data[[censored_variable]] == data[[censored_variable]][dim(data)[1]],censoring_indicator] <- 1
  }
  return(data)
}

set_last_as_censored <- function(data, censored_variable, censoring_indicator){
  if (!last_is_censored(data, censoring_indicator)){
    data[dim(data)[1],censoring_indicator] <- 0
    # repeat for duplicate values from bootstrapping
    data[data[[censored_variable]] == data[[censored_variable]][dim(data)[1]],censoring_indicator] <- 0
  }
  return(data)
}



