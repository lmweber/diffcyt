# Complete case analysis
#
# @param data 'data.frame'.
# @param censored_variable name of column containing censored data.
# @param censoring_indicator name of column containing indication if observed
#   (1) or censored (0) value in column 'censored_variable'.
# @param formula the formula for fitting the regression model (e.g. formula(y~x))
# @param regression_type The regression type to be used, one of 'lm' (linear
# regression), 'glm' (generalized linear regression), 'glmer' (generalized
# linear mixed-effects models). Default: 'lm'.
# @param weights name of column containing weights to be used in fitting the
#  regression model. Default = 'weights'. Ignored if no weights used.
# @param family The family to be used in the regression model.
#  Default = "binomial". Omitted if linear model ('lm') is used.
#
# @return list of elements 'data', 'betasMean' (mean regression coef of censored
#  covariate), 'betasVar' (mean variance of regression coef of censored covariate),
#  'fits' (regression fit)
# @export
#
# @examples
#  lm_formula <- formula(Y ~ Surv(X,I) + Z)
#  data <- simulate_singlecluster(100, lm_formula, type = "lm", n_levels_fixeff=2)
#  cc_out <- complete_case(data,"X","I",Y~X+Z)
#  summary(cc_out$fits)
complete_case <- function(data,
                          censored_variable,
                          censoring_indicator,
                          formula,
                          regression_type = c("lm", "glm", "glmer"),
                          weights = "weights",
                          family = "binomial"){
  regression_type <- match.arg(regression_type)
  # checking for right format of data, some type conversions
  data <- data_processing_for_imputation(data = data,
                                         censored_variable = censored_variable,
                                         censoring_indicator = censoring_indicator)
  
  # testing if enough observed events are present
  if (sum(data[[censoring_indicator]]) <= 3){
    stop("too few observed events for fitting regression models", call. = FALSE)
  }
  
  # create arguments list, only keep observed values
  args <- list(formula = formula,
               data = data[c(data[[censoring_indicator]] == 1), ])
  if (!identical(regression_type,"lm")) {
    args[["family"]] <- family
    if (!is.null(data[[weights]])) {
      args[["weights"]] <- data[[weights]][data[[censoring_indicator]] == 1]
    }
  }
  
  # regression fitting
  csi_fit <- do.call(regression_type, args)
  
  betas <- get_betas_from_fit(csi_fit, regression_type)
  vars <- get_variance_of_betas_from_fit(csi_fit)
  data_cmi <-
    list(
      data = data,
      betasMean = betas,
      betasVar = vars,
      fits = csi_fit
    )
  return(data_cmi)
}
