#' @importFrom stats runif 
kaplan_meier_imputation <- function(data, censored_variable, censoring_indicator, covariates = NULL){
  n <- dim(data)[1]
  censored_indices <- which(data[ ,censoring_indicator] == 0)
  # risk set calculation the same as in risk_set_imputation
  if (is.null(covariates)){
    Risk_Set <- purrr::map(censored_indices, ~ (.x+1):n)
  } else{
    Risk_Set <- risk_set_cov_adjusted(data, censored_variable, censoring_indicator, covariates)
  }
  # loop through each censored value
  est <- unlist(purrr::map_dbl(Risk_Set, function(j){
    subdata <- data[j,c(censored_variable,censoring_indicator)]
    subdata_sorted <- subdata[order(subdata[[censored_variable]]), ]
    # if last value is not observed set it to observed
    subdata_sorted <- set_last_as_observed(subdata_sorted, censored_variable, censoring_indicator)

    # fit kapplan meier curve
    km.fit <- summary(survival::survfit(survival::Surv(subdata_sorted[[censored_variable]], subdata_sorted[[censoring_indicator]]) ~ 1))
    # random draw according to survival time
    est_ind <- which(((1-km.fit$surv) >= runif(1)) == TRUE)[1]
    return(km.fit$time[est_ind])
  }))
  return(est)
}
