# impute censored values by adding a column to the data containing the completed
# data
# 
#' @importFrom magrittr %>%
#' @importFrom rlang :=
estimate_cens_vars <- function(data, censored_variable, censoring_indicator,
                              covariates = NULL, id = NULL,
                              imputation_method = c("mrl","rs","km","km_exp","km_wei","km_os"),
                              mi_reps = 1){
  imputation_method <- match.arg(imputation_method)
  n <- dim(data)[1]
  est_name <- paste0(censored_variable,"_est")
  # create an id column if missing
  if (is.null(id)){
    data$id <- seq_along(data[[censored_variable]])
    id <- "id"
  }

  # if last value is censored, set it to observed
  if (last_is_censored(data, censoring_indicator) & (imputation_method %in% c("mrl","rs","km"))) {
    last_censored <- TRUE
    data <- set_last_as_observed(data, censored_variable, censoring_indicator)
    # no censored values, return data
    if ( sum(data[[censoring_indicator]]) == length(data[[censoring_indicator]])){
      return(data)
    }
  } else{
    last_censored <- FALSE
  }
  
  # check that data is in correct format, some type conversions
  data_prep <- data_processing_for_imputation(data, censored_variable, censoring_indicator, covariates , id)

  # do estimation
  est <- switch(imputation_method,
                "mrl" = mean_residual_life_imputation(data_prep, censored_variable, censoring_indicator, covariates,id),
                "rs" = risk_set_imputation(data_prep, censored_variable, censoring_indicator, covariates,mi_reps),
                "km" = kaplan_meier_imputation(data_prep, censored_variable, censoring_indicator, covariates, mi_reps=mi_reps),
                "km_exp" =  kaplan_meier_imputation(data_prep, censored_variable, censoring_indicator, covariates, "exp",mi_reps),
                "km_wei" =  kaplan_meier_imputation(data_prep, censored_variable, censoring_indicator, covariates, "wei",mi_reps),
                "km_os" =  kaplan_meier_imputation(data_prep, censored_variable, censoring_indicator, covariates, "os",mi_reps)
  )

  if (mi_reps == 1) {
    duplicates_of_last_bool <- purrr::as_vector(data[n,id]) == purrr::as_vector(data[,id])
    
    # estimates of observed values are the same
    data[[est_name]] <- data[[censored_variable]]
    # estimates of censored values are real estimates
    data[c(data[[censoring_indicator]] == 0 & !duplicates_of_last_bool), est_name] <- est[, 1]
    
    # duplicates_of_last_bool <- FALSE
    # if last value was censored, set it back to censored
    if (last_censored){
      data <- set_last_as_censored(data, censored_variable, censoring_indicator)
      
      # for duplicates
      if (sum(duplicates_of_last_bool) > 1){
        data[n,est_name] <- purrr::as_vector(data[duplicates_of_last_bool,est_name])[1]
      }
    }
    return(data)  
  } else{
    return(asplit(est,2))
  }
}


# Conditional single imputation
#
# Imputes censored values according to \href{https://www.researchgate.net/publication/319246304_Improved_conditional_imputation_for_linear_regression_with_a_randomly_censored_predictor}{Atem et al. 2017}
#
# @param data 'data.frame'. Columns are the variables and
#  rows the samples.
# @param censored_variable name of column containing censored data
# @param censoring_indicator name of column containing indication if observed
#   (1) or censored (0) value in column 'censored_variable'
# @param response name of column containing the response (can be 'NULL')
# @param covariates name(s) of column(s) containing the covariate that influences
#  censoring
# @param id Default = 'NULL'. name of column containing id of sample. 
# @param imputation_method one of 'km','rs','mrl'. See \code{\link{conditional_multiple_imputation}}
# @param verbose Logical
# 
# @examples 
# lm_formula <- formula(Y ~ Surv(X,I) + Z)
# data <- simulate_singlecluster(50, lm_formula, type = "lm")
# conditional_single_imputation(
#   data = data,
#   censored_variable = "X",
#   censoring_indicator = "I",
#   response = "Y", 
#   covariates = "Z",
#   imputation_method = "km")
#' @importFrom magrittr %>%
conditional_single_imputation <- function(data, censored_variable,
                                          censoring_indicator, response = NULL,
                                          covariates = NULL, id = NULL,
                                          imputation_method = c("km","rs","mrl"),
                                          verbose = FALSE) {
  if (is.numeric(censored_variable)) censored_variable <- colnames(data)[[censored_variable]]
  if (is.numeric(censoring_indicator)) censoring_indicator <- colnames(data)[[censoring_indicator]]
  if (is.numeric(covariates) & !is.null(covariates)) covariates <- colnames(data)[covariates]
  if (is.null(id)){
    data$id <- seq_along(data[[censored_variable]])
    id <- "id"
  }
  if (is.numeric(id) & !is.null(id)) id <- colnames(data)[[id]]
  if (is.numeric(response) & !is.null(response)) response <- colnames(data)[[response]]
  imputation_method <- match.arg(imputation_method)
  # # check that data is in correct format, some type conversions, but keep covariates = NULL to not convert to factors
  data <- data_processing_for_imputation(data, censored_variable, censoring_indicator, covariates = NULL , id)
  
  est_name <- paste0(censored_variable,"_est")
  # check if censored values are present, if not return input
  nr_observed <- sum(data[[censoring_indicator]])
  length_data <- length(data[[censoring_indicator]])
  if ((nr_observed == length_data) |
      ((nr_observed == (length_data-1)) &
       # if last is censored cannot do calculations since it is transformed to 1
       (data[length_data,censoring_indicator] == 0)
      )
     ) {
    if (verbose) warning("No censored values, return input")
    data <- dplyr::mutate(data, !!est_name := data[[censored_variable]]) %>%
        dplyr::arrange(!!dplyr::sym(id)) 
    return(data)
  }
  # check if at least two observed values present
  if (sum(data[[censoring_indicator]]) < 2) {
    if (verbose) warning("Not enough observed values, return input")
    data <- dplyr::mutate(data, !!est_name := data[[censored_variable]]) %>%
      dplyr::arrange(!!dplyr::sym(id))
    return(data)
  }
  # estimates with or without cov,
  data <- estimate_cens_vars(data = data, censored_variable = censored_variable,censoring_indicator=censoring_indicator,
                      covariates = covariates, imputation_method = imputation_method, id = id)

  data <- dplyr::arrange(data, !!dplyr::sym(id)) 
  return(data)
}



