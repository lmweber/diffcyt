# fit cox proportional hazards model and return coefficients
coxph_coef_calc <- function(formula, data){
  cox_fit <- survival::coxph(formula, data = data)
  return(matrix(cox_fit$coefficients, ncol = 1))
}

# matrix multiplication between data and coxph coefficients with subsequent 
# centering and normalization
RS_cent_calc <- function(data, cox_coef_mat){
  RS <- as.matrix(data) %*%  cox_coef_mat
  return((RS-mean(RS))/stats::sd(RS))
}

# data needs to be sorted according to censored_variable
risk_set_cov_adjusted <- function(data, censored_variable, censoring_indicator,
                                  covariates, nn=10, wc=0.5){
  # check that correct weight is given
  if (wc < 0 | wc > 1){
    stop("weight 'wc' in 'risk_set_cov_adjusted' not between 0 and 1", 
         call. = FALSE)
  }
  
  # check that correct number of neighbors is given
  if (nn < 1){
    stop("number of neighbors to use 'nn' in 'risk_set_cov_adjusted' too low", 
         call. = FALSE)
  }
  
  n <- dim(data)[1]
  tmpformula <- as.formula(
    paste0("survival::Surv(data$", censored_variable, ", data$",
           censoring_indicator, ") ~ ",
           paste0(purrr::map_chr(covariates, ~ paste0("data$", .x)), collapse = "+")
    )
  )
  # make sure covariates are numeric
  data <- covariates_from_factor_to_numeric(data, covariates)
  censored_indices <- which(data[ ,censoring_indicator] == 0)
  event_indices <- which(data[ ,censoring_indicator] == 1)

  # cox proportional hazards coefficients
  cox_coef_c <- coxph_coef_calc(tmpformula, data[censored_indices, ])
  cox_coef_e <- coxph_coef_calc(tmpformula, data[event_indices, ])
  
  # if any coefficient is exactly zero return normal risk set (cannot calculate)
  if (anyNA(cox_coef_c) | anyNA(cox_coef_e) | any(cox_coef_c == 0) | any(cox_coef_e == 0)){ 
    return(purrr::map(censored_indices, ~ (.x+1):n))
  }
  
  # normalized, centered matrix
  RSc_cent <- RS_cent_calc(data[,covariates], cox_coef_c)
  RSe_cent <- RS_cent_calc(data[,covariates], cox_coef_e)
  
  # main distance calculation
  d2 <- sqrt(     wc * (outer(RSc_cent[,1], RSc_cent[,1], FUN = "-") ^ 2) +
              (1-wc) * (outer(RSe_cent[,1], RSe_cent[,1], FUN = "-") ^ 2)
            )
  rownames(d2) <- seq_len(n)
  colnames(d2) <- seq_len(n)
  
  # create risk set for all censored values
  Risk_Set <- purrr::map(censored_indices, function(j){
    startind <- which((data[[censored_variable]][j] < data[[censored_variable]]))[1]
    startind <- ifelse(is.na(startind), n, startind)
    if (n-startind >= 1){
      risk_set_j <- as.integer(names(sort(d2[j,startind:n])[seq_len(min(nn,n-startind+1))]))
    } else{
      risk_set_j <- n
    }
    return(risk_set_j)
  })
  return(Risk_Set)
}


risk_set_imputation <- function(data, censored_variable, censoring_indicator, covariates = NULL, mi_reps = 1){
  n <- dim(data)[1]
  # assumes sorted dataset according to censored variable
  censored_indices <- which(data[ ,censoring_indicator] == 0)
  if (is.null(covariates)){
    # risk set without covariates in sorted data
    Risk_Set <- purrr::map(censored_indices, ~ (.x+1):n)
  } else{
    # risk set with covariates in sorted data, see Hsu et al. 2006
    Risk_Set <- risk_set_cov_adjusted(data, censored_variable, censoring_indicator, covariates)
  }
  # estimate is a random value from the risk set
  est <-
    tryCatch(
      purrr::reduce(purrr::map(Risk_Set, ~ matrix(
        as.double(data[[censored_variable]][sample(.x, size = mi_reps, replace = TRUE)]), nrow = 1
      )), rbind),
      error = function(e)
        integer(0)
    )
  return(est)
}
