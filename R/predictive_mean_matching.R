# Predictive mean matching analysis
#
# Perform predictive mean matching by treating censored values as missing
#
# @param data 'data.frame'.
# @param censored_variable name of column containing censored data.
# @param censoring_indicator name of column containing indication if observed
#   (1) or censored (0) value in column 'censored_variable'.
# @param variables_for_imputation name of variable to use for imputation.
# @param formula the formula for fitting the regression model (e.g. formula(y~x)).
# @param regression_type The regression type to be used, one of 'lm' (linear
# regression), 'glm' (generalized linear regression), 'glmer' (generalized
# linear mixed-effects models). Default: 'lm'.
# @param repetitions number of repetitions for multiple imputation. Default: 10
# @param weights name of column containing weights to be used in fitting the
#  regression model. Default = 'weights'. Ignored if no weights used.
# @param family The family to be used in the regression model.
#  Default = "binomial". Omitted if linear model ('lm') is used.
#
# @return list of elements 'data', 'betasMean' (mean regression coef of censored
#  covariate), 'betasVar' (mean variance of regression coef of censored covariate),
#  'fits' (regression fits)
# @export
# @examples
#  lm_formula <- formula(Y ~ Surv(X,I) + Z)
#  data <- simulate_data(100, lm_formula, type = "lm", n_levels_fixeff=2)
#  pmm_out <- predictive_mean_matching(data,"X","I",c("Y","Z"),Y~X+Z)
#  summary(mice::pool(pmm_out$fits))
#' @importFrom utils capture.output
predictive_mean_matching <- function(data,
                                     censored_variable,
                                     censoring_indicator,
                                     variables_for_imputation = NULL,
                                     formula = NULL,
                                     regression_type = c("lm", "glm", "glmer"),
                                     repetitions = 10,
                                     weights = "weights",
                                     family = NULL){
  regression_type <- match.arg(regression_type)
  # checking for right format of data, some type conversions
  data <- data_processing_for_imputation(data = data,
                                   censored_variable = censored_variable,
                                   censoring_indicator = censoring_indicator)
  
  data_proc <- data
  # set censored values as missing
  data_proc[c(data_proc[[censoring_indicator]] == 0), censored_variable] <- NA
  
  # remove weights for pmm (linearly dependent on other covariate)
  if (!is.null(data[[weights]])){
    data_proc[[weights]] <- NULL
  }
  # set which variables to use for pmm, default is all
  if (!is.null(variables_for_imputation)){
    data_proc <- data_proc[ ,c(censored_variable, censoring_indicator, variables_for_imputation)]
  }
  # run mice, generating complete datasets
  capout <- capture.output({suppressWarnings(
    ini2 <- mice::mice(data = data_proc, m=repetitions, method = "pmm", maxit = 5)
  )})
  
  # construct arguments list
  args <- list(formula = formula)
  if (!identical(regression_type,"lm")) {
    args[["family"]] <- family
    if (!is.null(data[[weights]])) {
      args[["weights"]] <- data[[weights]]
    }
  }
  args_data <- mice::complete(ini2, action = "all")
  args_data <-  purrr::map(args_data, function(i){
    args[["data"]] <- i
    return(args)
  })
  
  # do fitting
  fits <- purrr::map(seq_along(args_data), ~ do.call(regression_type, args_data[[.x]]))
  
  data_cmi <- list(data = data,
                   betasMean = colMeans(get_betas_from_multiple_fits(fits,regression_type)),
                   betasVar = colMeans(get_variance_of_betas_from_multiple_fits(fits)),
                   fits = fits)
  return(data_cmi)
}
