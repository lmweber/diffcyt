# tail approx methods:
#   "lao" (last as observed): treat last values as being censored
#   "exp" (exponential): Brown-Hollander-Korwar completion of the product-limit 
#                       estimator, approximate tail by exponential distribution
#   "os" (order statistics): calculate expected order statistics from a fitted 
#                           weibull distribution. Replacement is then a random
#                           draw from one of those.
#   "wei" : fit weibull distribution, extrapolate
# 
# See: 'A Comparison of Several Methods of Estimating the Survival Function When 
#   There is Extreme Right Censoring' (M. L. Moeschberger and John P. Klein, 1985)
# 
#' @importFrom stats runif 
kaplan_meier_imputation <-
  function(data,
           censored_variable,
           censoring_indicator,
           covariates = NULL,
           tail_approx_method = c("lao", "exp", "os","wei"),
           mi_reps = 1) {
    tail_approx_method <- match.arg(tail_approx_method)
  n <- dim(data)[1]
  all_censored_bool <- data[ ,censoring_indicator] == 0
  # censored indices that have a lower value than the last observed one
  censored_indices <- which(all_censored_bool)# & 

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
  # censored indices that have a lower value than the last observed one
  censored_indices_with_risk_set <- censored_indices[!(censored_indices %in% censored_indices_no_risk_set)]
  
  # remove risk set of censored values that are higher than the last observed one (cannot fit survival curve with only censored data)
  if (length(censored_indices_no_risk_set) > 0){
    for (i in seq_along(censored_indices_no_risk_set)) {
      Risk_Set[[length(Risk_Set)]] <- NULL
    }
  }
  # additional requirements for tail approximation methods
  if (tail_approx_method == "exp"){
    # calculation of rate of exponential distribution
    km_fit_compl <- survival::survfit(survival::Surv(data[[censored_variable]], data[[censoring_indicator]]) ~ 1)
    km_fit_compl_summary <- summary(km_fit_compl)
    theta_compl <- -max(km_fit_compl_summary$time)/log(min(km_fit_compl_summary$surv))
  } else if (tail_approx_method == "os"){
    pot_replace_vals <- order_stat_surv_approx(data, censored_variable, censoring_indicator)
  } else if (tail_approx_method == "wei"){
    weibull_params <- params_from_weibull_fit(data, censored_variable, censoring_indicator)
  } else if (tail_approx_method == "lao"){
    # if last value is not observed set it to observed
    data <- set_last_as_observed(data , censored_variable, censoring_indicator)
  }
  
  # loop through each censored value
  est <- purrr::map(seq_along(censored_indices_with_risk_set) , function(js){
    current_cens_val <- purrr::as_vector(data[censored_indices_with_risk_set[js],censored_variable])
    
    subdata_sorted <- data[sort(Risk_Set[[js]]),c(censored_variable,censoring_indicator)]
    
    # fit kapplan meier curve
    km_fit <- survival::survfit(survival::Surv(subdata_sorted[[censored_variable]], subdata_sorted[[censoring_indicator]]) ~ 1)
    km_fit_summary <- summary(km_fit)
    # random draw according to survival time
    replacement_value <- surv_tail_lao(km_fit_summary, mi_reps)
    # check if some values that are inferred are higher than the highest observed one
    too_high <- replacement_value >= km_fit_summary$time[which.min(km_fit_summary$surv)] &
      km_fit_summary$n.event[which.min(km_fit_summary$surv)] < 1
    # for each inferred value that is higher than the highest observed one apply special procedure
    if (sum(too_high) > 0){
      if(tail_approx_method == "exp"){
        # random draw from truncated exponential
        x <- runif(sum(too_high),min = pexp(current_cens_val, theta_compl))
        replacement_value[too_high] <- qexp(x, theta)
      } else if (tail_approx_method == "wei"){
        # random draw from truncated weibull
        x <- runif(sum(too_high),min = pweibull(current_cens_val, weibull_params$shape, weibull_params$scale))
        replacement_value[too_high] <- qweibull(x, weibull_params$shape, weibull_params$scale)
      }  else if (tail_approx_method == "os"){
        # random draw from order statistics
        pot_replace_vals_adj <- pot_replace_vals[pot_replace_vals >= current_cens_val]
        replacement_value[too_high] <- matrix(resample(pot_replace_vals_adj,sum(too_high),replace = TRUE),nrow=1)
      }
    }
    return(replacement_value)
  }) %>% purrr::reduce(rbind)
  
  # for censored values higher than the highest observed one impute with exponential distribution
  if (length(censored_indices_no_risk_set) > 0){
      est <- rbind(est,purrr::map(censored_indices_no_risk_set, function(i){
        current_cens_val <- purrr::as_vector(data[i,censored_variable])
          if(tail_approx_method == "exp"){
            # random draw from truncated exponential
            x <- runif(mi_reps,min = pexp(current_cens_val, theta_compl))
            matrix(qexp(x, theta_compl),nrow = 1)
          } else if (tail_approx_method == "os"){
            pot_replace_vals_adj <- pot_replace_vals[pot_replace_vals >= current_cens_val]
            mat <- matrix(matrix(resample(pot_replace_vals_adj,mi_reps,replace = TRUE),nrow=1),nrow = 1)
            return(mat)
          } else if (tail_approx_method == "wei"){
            # random draw from truncated weibull
            x <- runif(mi_reps,min = pweibull(current_cens_val, weibull_params$shape, weibull_params$scale))
            matrix(qweibull(x, weibull_params$shape, weibull_params$scale),nrow = 1)
          } else if (tail_approx_method == "lao"){
            matrix(rep(data[max(censored_indices_no_risk_set),censored_variable],mi_reps),nrow=1)
          }
      }) %>% purrr::reduce(rbind))
    }
  dimnames(est) <- NULL
  return(est)
  }

# resample function from 'sample' manual to get correcto behaviour if length of x is 1
resample <- function(x, ...) x[sample.int(length(x), ...)]


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
    return(matrix(
    resample(
      km_fit_summary$time,
      size = mi_reps,
      replace = TRUE,
      prob = prob
    ),
    nrow = 1
  ))
  }
}




dtweibull <- function(x, shape, scale = 1, left_trunc=0) {
  shape/scale*(x/scale)^(shape-1)*exp(-((x/scale)^shape-(left_trunc/scale)^shape))
}

ptweibull <- function(x, shape, scale = 1, left_trunc=0) {
  1-exp(-((x/scale)^shape-(left_trunc/scale)^shape))
}

integrand <- function(x,r,n,shape=1, scale=1, left_trunc=0) {
  x * (1 - ptweibull(x, shape, scale,left_trunc))^(r-1) * ptweibull(x, shape, scale,left_trunc)^(n-r) * dtweibull(x, shape, scale,left_trunc)
}

expected_j_order_stat <- function(r,n, shape=1, scale=1, left_trunc=0) {
  (1/beta(n-r+1,n-(n-r+1)+1)) * integrate(integrand,left_trunc,Inf, n-r+1, n, shape, scale, left_trunc)$value
}

params_from_weibull_fit <- function(data, censored_variable, censoring_indicator){
  tmpdat <- data.frame(left=data[[censored_variable]], 
                       right = ifelse(data[[censoring_indicator]]==1,data[[censored_variable]],rep(NA,length(data[[censored_variable]]))))
  tmpdat <- tmpdat[order(tmpdat$left),]
  return(as.list(coef(fitdistrplus::fitdistcens(tmpdat,"weibull"))))
}

  
order_stat_surv_approx <- function(data, censored_variable, censoring_indicator){
  all_censored_bool <- data[ ,censoring_indicator] == 0
  # max observed value
  max_obs <- max(purrr::as_vector(data[data[[censoring_indicator]] == 1,censored_variable]))
  # censored indices that have a higher value than the last observed one
  censored_indices_no_risk_set <- which(all_censored_bool & 
                                          data[ ,censored_variable] > max_obs)
  # fit weibull to data to get parameters
  weibull_params <- params_from_weibull_fit(data, censored_variable, censoring_indicator)
  # calculate order statistics of censored values higher than highest observed one
  return(sapply(censored_indices_no_risk_set, function(i){
    expected_j_order_stat(i,dim(data)[1],weibull_params$shape,weibull_params$scale,left_trunc = max_obs)
  }))
}





