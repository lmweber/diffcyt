#' @importFrom stats runif 
kaplan_meier_imputation <- function(data, censored_variable, censoring_indicator, covariates = NULL, tail_approx_method = "lao", mi_reps = 1){
  n <- dim(data)[1]
  all_censored_bool <- data[ ,censoring_indicator] == 0
  # censored indices that have a lower value than the last observed one
  censored_indices <- which(all_censored_bool)# & 
                              # data[ ,censored_variable] <= max(purrr::as_vector(data[data[[censoring_indicator]] == 1,censored_variable])))
  
  km_fit_compl <- survival::survfit(survival::Surv(data[[censored_variable]], data[[censoring_indicator]]) ~ 1)
  km_fit_compl_summary <- summary(km_fit_compl)
  theta_compl <- -max(km_fit_compl_summary$time)/log(min(km_fit_compl_summary$surv))
    # risk set calculation the same as in risk_set_imputation
  if (is.null(covariates)){
    Risk_Set <- purrr::map(censored_indices, ~ (.x+1):n)
  } else{
    Risk_Set <- risk_set_cov_adjusted(data, censored_variable, censoring_indicator, covariates)
  }
  # max observed value
  max_obs <- max(purrr::as_vector(data[data[[censoring_indicator]] == 1,censored_variable]))
  # censored indices that have a higher value than the last observed one
  censored_indices_no_risk_set <- which(all_censored_bool & 
                                          data[ ,censored_variable] > max_obs)
  # remove risk set of censored values that are higher than the last observed one (cannot fit survival curve with only censored data)
  if (length(censored_indices_no_risk_set) > 0){
    for (i in seq_along(censored_indices_no_risk_set)) {
      Risk_Set[[length(Risk_Set)]] <- NULL
    }
  }
  # loop through each censored value
  est <- purrr::map(Risk_Set, function(js){
    subdata <- data[js,c(censored_variable,censoring_indicator)]
    subdata_sorted <- subdata[order(subdata[[censored_variable]]), ]
    
    if (tail_approx_method == "lao"){
      # if last value is not observed set it to observed
      subdata_sorted <- set_last_as_observed(subdata_sorted, censored_variable, censoring_indicator)
    }
    
    # fit kapplan meier curve
    km_fit <- survival::survfit(survival::Surv(subdata_sorted[[censored_variable]], subdata_sorted[[censoring_indicator]]) ~ 1)
    km_fit_summary <- summary(km_fit)
    # random draw according to survival time
    replacement_value <- switch (tail_approx_method,
                                 "lao" = surv_tail_lao(km_fit_summary,mi_reps),
                                 "exp" = surv_tail_exp(km_fit_summary,mi_reps,theta_compl)
    )
    return(replacement_value)
  }) %>% purrr::reduce(rbind)
  
  # for censored values higher than the highest observed one impute with exponential distribution
  if (length(censored_indices_no_risk_set) > 0){
    if (tail_approx_method == "exp"){
      est <- rbind(est,purrr::map(censored_indices_no_risk_set, function(x){
        # random draw from truncated exponential
        x <- runif(mi_reps,min = pexp(max(km_fit_compl_summary$time), theta_compl))
        matrix(qexp(x, theta_compl),nrow = 1)
      }) %>% purrr::reduce(rbind))
    }
  }
  
  return(est)
}

# random replacement from survival curve
surv_tail_lao <- function(km_fit_summary,mi_reps=1){
  if(length(km_fit_summary$surv) == 1){
    return(km_fit_summary$time)
  } else {
    prob <- diff(c(0,1 - km_fit_summary$surv, 1))
    prob[length(prob)-1] <- sum(prob[(length(prob)-1):length(prob)])
    prob <- prob[seq_len(length(prob)-1)]
    if(length(prob) == 0){
      prob <- 1
    }
    # print(km_fit_summary)
    # print(prob)
    return(matrix(
    sample(
      km_fit_summary$time,
      size = mi_reps,
      replace = TRUE,
      prob = prob
    ),
    nrow = 1
  ))
  }
}


surv_tail_exp <- function(km_fit_summary,mi_reps = 1, theta = -max(km_fit_summary$time)/log(min(km_fit_summary$surv))){
  # normal replacement if random draw from observed region
  replacement_value <- surv_tail_lao(km_fit_summary,mi_reps)
  # values that are higher than observed
  too_high <- replacement_value >= km_fit_summary$time[which.min(km_fit_summary$surv)] &
    km_fit_summary$n.event[which.min(km_fit_summary$surv)] < 1
  if (sum(too_high) > 0){
    # random draw from truncated exponential
    x <- runif(sum(too_high),min = pexp(max(km_fit_summary$time), theta))
    replacement_value[too_high] <- qexp(x, theta)
  }
  return(replacement_value)
}

