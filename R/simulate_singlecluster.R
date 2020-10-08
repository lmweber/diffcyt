#' Simulation of data with a censored covariate
#'
#' Function to simulate an association between a censored covariate and a predictor.
#'
#' @param n number of samples
#' 
#' @param formula the formula to specify the structure in the data. The censored
#'  variable should be written in the following format: 'Surv(X,I)', where 'X' is the observed
#'  value, and 'I' is the event indicator (1 if observed, 0 if censored).
#'  A full example is: 'Y ~ Surv(X,I) + Covariate + (1|Random_effect)'.
#'  
#' @param type which regression type is used, one of 'lm', 'glm', 'glmer'. For
#'  the generalized linear models the response is binomial with a logistic link
#'  function. default = 'lm'.
#'  
#' @param b the regression coefficients, either
#' \describe{
#'  \item{NUll}{will us 0 for the intercept and 1
#'  for the remaining coefficients} 
#'  \item{a vector with regression coefficients}{the length has to be (1 
#'  (intercept) + number of covariates (including the censored covariate))}
#'   }
#'   
#' @param n_levels_fixeff The number of levels to use for each covariate, e.g.
#'  for two covariates: c(10,100). If NULL sets all to 2 (two groups).
#'  
#' @param n_levels_raneff The number of levels to use for each random effect. 
#' If NULL sets to 'n' (observation level random effects).
#' 
#' @param weibull_params The parameters for the distribution of the censored
#' variable and the censoring time. Should be a list of lists, where the elements
#' of the outer lists are 'X' the true value and 'C' the censoring time. The
#' inner lists should have two keywords, 'shape' and 'scale', for the parameters
#' of the Weibull distribution (See \code{\link[stats]{rweibull}}).
#' 
#' @param censoring_dependent_on_covariate Logical. If
#'  censoring should depend on a covariate. The respective covariate needs to
#'  have only two levels ('n_level_fixeff'=2). Will use first covariate
#'  in formula.
#'  
#' @param weibull_params_covariate_dependent_censoring list with two elements,
#'  shape and scale, representing the parameters of a weibull distribution for
#'  the second level of a covariate if 'censoring_dependent_on_covariate'=TRUE.
#'  
#' @param error_variance positive double. Variance of additional gaussian
#'  noise to add in the linear sum of the predictors. For linear regression
#'  this is the only error added. Otherwise it should be set to zero. default = 0.
#'  
#' @param variance_raneff positive double vector of the length of 
#'  'n_levels_raneff'. The variance of the gaussian distributed
#'  random effect covariates. default = 0.5.
#'  
#' @param transform_fn function to transform censored covariate or one of
#'  'identity' (no transformation), 'boxcox' (box-cox transformation),
#'  'boxcox_positive' (box-cox transformation and translation to all positive
#'   values), 'log_positive' (log transformation and translation to all positive
#'   values). The transformation is applied before the response is modeled.  
#'   default = 'identity'.
#'   
#' @param verbose verbose
#' 
#' @export
#' 
#' @importFrom stats runif rnorm rweibull
#' @importFrom rlang :=
#' 
#' @return \code{\link[tibble]{tibble}} 
#' @examples
#' # single differential cluster
#'  glmer_formula <- formula(Y ~ Surv(X,I) + Z + (1|R))
#'  simulate_singlecluster(100, glmer_formula, type = "glmer")
#'  
simulate_singlecluster <- function(n,
                          formula,
                          type = c("lm","glm","glmer"),
                          b = NULL,
                          n_levels_fixeff = NULL,
                          n_levels_raneff = NULL,
                          weibull_params = list(X = list(shape = 0.5, scale = 0.25),
                                                C = list(shape = 1, scale = 0.25)),
                          censoring_dependent_on_covariate = FALSE,
                          weibull_params_covariate_dependent_censoring = list(shape = 1, scale = 0.1),
                          error_variance = 0,
                          variance_raneff = 0.5,
                          transform_fn = "identity",
                          verbose = FALSE){
  ######  input checks
  
  # sample size
  if ( (n < 1) | ((n %% 1) != 0)) stop("not valid number of samples 'n'",call. = FALSE)
  
  # formula
  uncensored_formula <- create_glmm_formula(formula)
  if (verbose) {
    formula_as_character <- as.character(uncensored_formula)
    message(cat("uncensored Formula: ",formula_as_character[2],formula_as_character[1],formula_as_character[3],"\n"))
  }
  type <- match.arg(type)
  
  # variables from formula
  variables <- extract_variables_from_formula(formula)
  # if (is.null(number_of_differential_clusters)){
    number_of_differential_clusters <- 1
  # }
  if (verbose) message(cat("Formula elements: \n",paste("\t",names(variables),"\t", variables, collapse = "\n"),"\n"))
  
  # regression coefficients
  if (is.null(b)){
    b <-
      rep(list(b = c(0, rep(
        1, 1 + length(variables$covariates)
      ))), times = number_of_differential_clusters)
  } else if (is.list(b) & length(b)==1){
    b <- rep(b,number_of_differential_clusters)
  } else if (!is.list(b) & length(b)==(2 + length(variables$covariates))){
    b <- rep(list(b=b),number_of_differential_clusters)
  }
  stopifnot(is.list(b) & (length(b) == number_of_differential_clusters))
  
  # weibull parameters
  stopifnot(all(unlist(weibull_params) >= 0))
  stopifnot(all(unlist(weibull_params_covariate_dependent_censoring) >= 0))
  
  # covariate dependent censoring needs two levels of the covariate
  if(censoring_dependent_on_covariate) { 
    stopifnot(n_levels_fixeff[1] == 2 | is.null(n_levels_fixeff))
  }
  if (!is.null(variables$covariates) & is.null(n_levels_fixeff)){
    n_levels_fixeff <- rep(2, length(variables$covariates))
  }
  if (!is.null(variables$covariates) & length(n_levels_fixeff)==1){
    n_levels_fixeff <- rep(n_levels_fixeff, length(variables$covariates))
  }
  if (!is.null(variables$random_covariates) & is.null(n_levels_raneff)){
    n_levels_raneff <- rep(n, length(variables$random_covariates))
  }
  stopifnot((length(n_levels_fixeff) == length(variables$covariates)))
  stopifnot((length(n_levels_raneff) == length(variables$random_covariates)))
  

  # fixed effects covariates
  if (!is.null(variables$covariates)){ 
    fixed_effects <- lapply(seq_along(n_levels_fixeff), function(x) gl(n = n_levels_fixeff[x], k = n/n_levels_fixeff[x]))
  }
  
  # need random covariates if type is glmer
  stopifnot(!is.null(variables$random_covariates) | (type != "glmer")) 
  
  # random effects
  if (!is.null(variables$random_covariates) & (type == "glmer")){ 
    if (length(variance_raneff) != length(variables$random_covariates) & length(variance_raneff) == 1){
      variances_random_intercepts <- rep(variance_raneff, length(n_levels_raneff))
    } else {
      variances_random_intercepts <- variance_raneff
    }
  }
  
  
  ###### start data simulation part

  # True survival times
  X1 <-  rweibull(n, weibull_params$X$shape, weibull_params$X$scale) 
  # potential censoring times
  C <- rweibull(n, weibull_params$C$shape, weibull_params$C$scale) 
  # if censoring depends on covariate 
  if (!is.null(variables$covariates) & censoring_dependent_on_covariate){
    C2 <- rweibull(n, weibull_params_covariate_dependent_censoring$shape,
                   weibull_params_covariate_dependent_censoring$scale)
    if ((censoring_dependent_on_covariate %in% variables$covariates) &
        !is.logical(censoring_dependent_on_covariate)){
      covariate_for_censoring <- which(censoring_dependent_on_covariate == variables$covariates)
    } else {
      covariate_for_censoring <- 1
    }
    C <- ifelse(levels(fixed_effects[[covariate_for_censoring]])[1] == fixed_effects[[covariate_for_censoring]], C, C2)
  }
  
  # observed survival times
  O1 <- pmin(X1,C) 
  # indicator vector, 1 if event happened
  I1 <- (X1 < C) * 1 
  if (verbose) cat("Censored values: ", length(I1)-sum(I1), " / ", length(I1),"\n")
  # transform X1 and O1
  X1_O1 <- do_transformation(c(X1,O1),transform_fn,verbose)
  X1 <- X1_O1[seq_along(X1)]
  O1 <- X1_O1[seq_along(X1)+length(X1)]

  # samples sizes (nr. cells) per patient
  size_tot <- round(runif(n,1e4, 1e5)) 
  
  # put data together
  col_data <- dplyr::tibble(ID = seq_len(n),
                       TrVal = X1,
                       !!variables$censored_variable := O1,
                       !!variables$censoring_indicator := I1,
                       n_cells = size_tot)

  # covariates list, rep(1,..) is for later matrix multiplication of intercept
  covariates_list <- list(rep(1,length(X1)),X1)
  
  # calculate actual fixed effects covariates
  if (!is.null(variables$covariates)){ 
    values_fixed_effects <- lapply(seq_along(variables$covariates), function(x){
      # random value for each sample with the same covariate level
      pat_intercept <- c(0, seq_len(n_levels_fixeff[x]-1)) 
      rep(pat_intercept, each = n/n_levels_fixeff[x])
    })
    covariates_list <- append(covariates_list,lapply(seq_along(variables$covariates), function(x){values_fixed_effects[[x]]}))
    
    # add to output data
    for (i in seq_along(variables$covariates)) {
      col_data <- col_data %>% dplyr::mutate(!!variables$covariates[i] := values_fixed_effects[[i]])
    }
  }
  
  # random effects covariates
  if (!is.null(variables$random_covariates) & (type == "glmer")){ 
    values_random_intercepts <- lapply(seq_along(n_levels_raneff), function(x){
      # random intercept for each sample with the same covariate level
      pat_intercept <- rnorm(n_levels_raneff[x], 0, variances_random_intercepts[x]) 
      rep(pat_intercept, each = n/n_levels_raneff[x])
    })
    list_random_covariates <- lapply(seq_along(variables$random_covariates), function(x){values_random_intercepts[[x]]})
    
    # add to output data
    for (i in seq_along(variables$random_covariates)) {
      col_data <- col_data %>% dplyr::mutate(!!variables$random_covariates[i] := values_random_intercepts[[i]])
    }
  } else{
    list_random_covariates <- NULL
  }
  
  # create 'number_of_differential_clusters' data sets
  out_counts_ls <- purrr::map(seq_len(number_of_differential_clusters), function(i){
    # calculate linear sum with covariates, regression coefficients and random effects
    linear_sum <- linear_sum(covariates_list = covariates_list,
                             betas = b[[i]],
                             random_covariates_list = list_random_covariates)
    
    # linear regression
    if (type == "lm") {
      # add random error in linear regression
      Y <-  linear_sum + rnorm(n,0,error_variance) 
      # True dependency
      Y_True <- linear_sum
      
    # logistic regression
    }else if (type %in% c("glm","glmer")) { 
      # True dependency
      Y_True <- (1/(1+exp(-linear_sum)))
      # possibility to set additional normal error in logistic regression
      Y_with_error <- Y_True + rnorm(length(Y_True),0,error_variance)
      Y_with_error <- ifelse(Y_with_error>1 | Y_with_error < 0, Y_True, Y_with_error)
      # binomial link function
      Y_size <-  unlist(purrr::map2(Y_with_error,size_tot, ~ rbinom(1, .y, .x )))
      # proportion, further used
      Y <- Y_size/size_tot
    }
    
    # response name to use as label
    tmp_response_name <- variables$response
    
    y_tib <- tibble::tibble(!!tmp_response_name := as.vector(Y),
                            !!paste0(tmp_response_name,"_True") := as.vector(Y_True))
    
    return(list(col_data = y_tib))
  })
  
  
  # for single cluster data
  return(dplyr::bind_cols(col_data,out_counts_ls[[1]]$col_data))
}
