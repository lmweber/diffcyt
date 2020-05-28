# mean residual life imputation, aka conditional imputation
# 
#' @importFrom magrittr %>%
#' @importFrom rlang :=
mean_residual_life_imputation <- function(data, 
                                          censored_variable, 
                                          censoring_indicator, 
                                          covariates = NULL,
                                          id = NULL){
  n <- dim(data)[1]
  if (is.null(id)){
    data$id <- seq_along(data[[censored_variable]])
    id <- "id"
  }
  
  # conditional single imputation from Atem et al. 2017 without covariates
  if (is.null(covariates)){
    # fit kapplan meier survival curve
    km_fit <- survival::survfit(survival::Surv(data[[censored_variable]], data[[censoring_indicator]]) ~ 1)
    data$S <- summary(km_fit, times = data[[censored_variable]])$surv
    
    # only censored values
    tmp_censored <- purrr::as_vector(data[c(data[[censoring_indicator]] == 0), censored_variable])
    
    # create matrix with two times the same censored variable column but shifted
    # by one. for more efficient calculation later on
    data_matrix_cens_var <- as.matrix(cbind(data[-n, censored_variable],data[-1, censored_variable]))
    # the same with the survival times
    data_matrix_surv <- as.matrix(cbind(data[-n, "S"],data[-1, "S"]))
    
    # calculation of numerator from equation
    numerator <- purrr::map_dbl(tmp_censored, function(tmp_cens_val) {
      tmpinds <- which(data_matrix_cens_var[, 1] > tmp_cens_val)
      if (!(length(tmpinds) == 1)){
        sum((rowSums(data_matrix_surv[tmpinds, ]) *
               (data_matrix_cens_var[tmpinds, 2] - data_matrix_cens_var[tmpinds, 1])))
      }else{
        sum((sum(data_matrix_surv[tmpinds, ]) *
               diff(data_matrix_cens_var[tmpinds, ])))
      }
    })
    
    # assemble another matrix
    tmp_E0 <- data[c(data[[censoring_indicator]] == 0), c("S",censored_variable)]
    tmp_E0$x <- numerator
    tmp_E0 <- as.matrix(tmp_E0)
    
    # final calculation, get estimates for all censored values
    est <- tmp_E0[ ,3] / (2 * tmp_E0[ ,1]) + tmp_E0[ ,2]

    
  # conditional single imputation from Atem et al. 2017 with covariates
  } else{
    # formula for cox prop haz model
    tmpformula <- as.formula(
      paste0("survival::Surv(data$", censored_variable, ", data$",
             censoring_indicator, ") ~ ",
             paste0(purrr::map_chr(covariates, ~ paste0("data$", .x)), collapse = "+")
      )
    )
    # make sure covariates are numeric
    data <- covariates_from_factor_to_numeric(data, covariates)
    
    # fit cox proportional hazards model
    cox_fit <- survival::coxph(tmpformula, data = data)
    cox_coef_z <- cox_fit$coefficients
    data <- suppressMessages(
      dplyr::inner_join(data, 
                        dplyr::as_tibble(survival::basehaz(cox_fit)) %>% 
                          dplyr::rename(!!censored_variable := "time")))
    data <- data %>% dplyr::mutate(S = exp(-hazard)) %>% dplyr::select(-hazard)

    # data subset only censored entries
    id_cens <- as.matrix(data[c(data[ , censoring_indicator] == 0), c(id, covariates,censored_variable)])

    # for matrix calculations
    data_matrix <- as.matrix(data[, c("S",censored_variable)])
    # for matrix calculation of data_matrix[x,]-data_matrix[x-1,]
    data_matrix_cens_var <- as.matrix(cbind(data[-n, censored_variable],data[-1, censored_variable]))
    # loop through each censored entry, calculation of numerator of equation
    numerator <- purrr::map_dbl(id_cens[ ,id], function(C){
      # exponent calculation
      data_matrix_exponent <- data_matrix[, 1] ^ exp(sum(cox_coef_z * id_cens[which(id_cens[ , id] == C) , covariates]))
      data_matrix_exponent <- cbind(data_matrix_exponent[-n],data_matrix_exponent[-1])
      # all values with higher censored variable values than the current one
      which_cens <- which(data_matrix_cens_var[, 1] > as.double(id_cens[which(id_cens[ , id] == C)[1], censored_variable]))
      # normal case
      if (!(length(which_cens) == 1)){
        sum((rowSums(data_matrix_exponent[which_cens, ]) *
               (data_matrix_cens_var[which_cens, 2] - data_matrix_cens_var[which_cens, 1])))
        # different calculation needed if only one value present
      }else{
        sum((sum(data_matrix_exponent[which_cens, ]) *
               diff(data_matrix_cens_var[which_cens, ])))
      }
    })
    # combine to matrix
    tmp_E0 <- cbind(data[c(data[[censoring_indicator]] == 0), c("S",censored_variable,id)],numerator)
    tmp_E0 <- as.matrix(tmp_E0)
    tmpind <- unlist(lapply(seq_along(tmp_E0[ ,3]) , function(x) {which(tmp_E0[x,3] == data[ , id])[1]}))
    tmp_cov_mat <- as.matrix(data[tmpind, covariates])
    
    # end calculation, get estimates for all censored values
    est <- tmp_E0[,4] / (2 * tmp_E0[,1] ^ exp(rowSums(t(cox_coef_z*t(tmp_cov_mat))))) +  tmp_E0[ ,2]
  }
  return(matrix(est,ncol=1))
}
