#' data simulation function
#'
#' Simulation of data with a censored covariate
#'
#' @param n number of samples
#' @param formula the formula to specify the structure in the data. The censored
#'  variable should be written as following: 'Surv(X,I)', where 'X' is the observed
#'  value, and 'I' is the censoring indicator (1 if observed, 0 if censored).
#'  A full example is: 'Y ~ Surv(X,I) + Covariate + (1|Random_effect)'.
#' @param type which regression type is used, one of 'lm', 'glm', 'glmer'. For
#'  the generalized linear models the response is binomial with a logistic link
#'  function
#' @param b the regression coefficients, if Null uses 0 for the intercept and 1
#'  for the remaining coefficients
#' @param n_levels_fixeff The number of levels to use for each covariate, e.g.
#'  for two covariates: c(10,100). If null uses 'n' meaning all samples differ
#'  from each other (if random effect meaning observation level).
#' @param n_levels_raneff The number of levels to use for each random effect
#' @param weibull_params The parameters for the distribution of the censored
#' variable and the censoring time. Should be a list of lists, where the elements
#' of the outer lists are 'X' the true value and 'C' the censoring time. The
#' inner lists should have two keywords, 'shape' and 'scale', for the parameters
#' of the Weibull distribution.
#' @param number_of_clusters Positive Integer. The number of clusters per true differential
#'  cluster for testing \code{\link{testDA_censoredGLMM}}. The total number of clusters is
#'  'number_of_clusters' * 'number_of_differential_clusters'. If NULL (default)
#'  only one cluster is used (running \code{\link{conditional_multiple_imputation}})
#' @param number_of_differential_clusters Positive Integer. Total number of clusters
#'  with a true signal.
#' @param censoring_dependent_on_covariate Logical or name of covariate. If
#'  censoring should depend on a covariate. The respective covariate needs to
#'  have only two levels ('n_level_fixeff'=2).
#' @param weibull_params_covariate_dependent_censoring list with two elements,
#'  shape and scale, representing the parameters of a weibull distribution for
#'  the second level of a covariate if 'censoring_dependent_on_covariate'=TRUE.
#' @param error_variance positive double. Variance of additional gaussian
#'  noise to add in the linear sum of the predictors. For linear regression
#'  this is the only error added.
#' @param variance_fixeff positive double. The variance of the gaussian distributed
#'  fixed effect covariates
#' @param variance_raneff positive double. The variance of the gaussian distributed
#'  random effect covariates
#' @param transform_fn function to transform censored covariate or one of
#'  'identity' (no transformation), 'boxcox' (box-cox transformation),
#'  'boxcox_positive' (box-cox transformation and translation to all positive
#'   values), 'log_positive' (log transformation and translation to all positive
#'   values). default = 'identity'
#' @param verbose verbose
#' @return simulated dataset or list of simulated dataset if
#'  'number_of_clusters' != NULL
#' @examples
#'  lm_formula <- formula(Y ~ Surv(X,I) + Z)
#'  simulate_data(100, lm_formula, type = "lm")
#'
#'  glm_formula <- formula(Y ~ Surv(X,I) + Z)
#'  simulate_data(100, glm_formula, type = "glm")
#'
#'  glmer_formula <- formula(Y ~ Surv(X,I) + Z + (1|R))
#'  simulate_data(100, glmer_formula, type = "glmer")
#' @export
simulate_data <- function(n,
                          formula,
                          type = c("lm","glm","glmer"),
                          b = NULL,
                          n_levels_fixeff = NULL,
                          n_levels_raneff = NULL,
                          weibull_params = list(X = list(shape = 0.75, scale = 2),
                                                C = list(shape = 0.5, scale = 10)),
                          number_of_clusters = NULL,
                          number_of_differential_clusters = NULL,
                          censoring_dependent_on_covariate = FALSE,
                          weibull_params_covariate_dependent_censoring = list(shape = 1, scale = 10),
                          error_variance = 0,
                          variance_fixeff = 0.5,
                          variance_raneff = 0.5,
                          transform_fn = "identity",
                          verbose = FALSE){
  type = match.arg(type)
  if (n<1) stop("negative length 'n'",call. = FALSE)
  uncensored_formula <- create_glmm_formula(formula)
  if (verbose) {
    formula_as_character <- as.character(uncensored_formula)
    message(cat("uncensored Formula: ",formula_as_character[2],formula_as_character[1],formula_as_character[3],"\n"))
  }
  variables <- extract_variables_from_formula(formula)
  if (is.null(number_of_differential_clusters)){
    number_of_differential_clusters <- 1
  }
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
  stopifnot(all(unlist(weibull_params) >= 0))
  stopifnot(all(unlist(weibull_params_covariate_dependent_censoring) >= 0))
  if(censoring_dependent_on_covariate) stopifnot(n_levels_fixeff == 2 | is.null(n_levels_fixeff))
  if (verbose) message(cat("Formula elements: \n",paste("\t",names(variables),"\t", variables, collapse = "\n"),"\n"))
  if (!is.null(variables$covariates) & is.null(n_levels_fixeff)){
    n_levels_fixeff <- rep(n, length(variables$covariates))
  }
  if (!is.null(variables$covariates) & length(n_levels_fixeff)==1){
    n_levels_fixeff <- rep(n_levels_fixeff, length(variables$covariates))
  }
  if (!is.null(variables$random_covariates) & is.null(n_levels_raneff)){
    n_levels_raneff <- rep(n, length(variables$random_covariates))
  }
  stopifnot((length(n_levels_fixeff) == length(variables$covariates)))
  stopifnot((length(n_levels_raneff) == length(variables$random_covariates)))
  stopifnot(!is.null(variables$random_covariates) | (type != "glmer")) # need random covariates if type is glmer
  if (!is.null(variables$covariates)){ # fixed effects covariates
    fixed_effects <- lapply(seq_along(n_levels_fixeff), function(x) gl(n = n_levels_fixeff[x], k = n/n_levels_fixeff[x]))
  }



  X1 <-  rweibull(n, weibull_params$X$shape, weibull_params$X$scale) # True survival times
  C <- rweibull(n, weibull_params$C$shape, weibull_params$C$scale) # potential censoring times
  if (!is.null(variables$covariates) & censoring_dependent_on_covariate){ # if censoring depends on covariate
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
  O1 <- pmin(X1,C) # observe survival times
  I1 <- (X1 < C) * 1 # indicator vector, 1 if event happened
  if (verbose) cat("Censored values: ", length(I1)-sum(I1), " / ", length(I1),"\n")
  # transform X1 and O1
  if(verbose) cat("\nTrVal:\t")
  X1 <- do_transformation(X1,transform_fn,verbose)
  if(verbose) cat("X:\t")
  O1 <- do_transformation(O1,transform_fn,verbose)

  size_tot <- round(runif(n,1e4, 1e5)) # samples sizes (nr. cells) per patient
  out <- dplyr::tibble(ID = seq_len(n),
                       TrVal = X1,
                       !!variables$censored_variable := O1,
                       !!variables$censoring_indicator := I1,
                       size_tot = size_tot)


  covariates_list <- list(rep(1,length(X1)),X1)
  if (!is.null(variables$covariates)){ # fixed effects covariates
    fixed_effects <- lapply(seq_along(n_levels_fixeff), function(x) gl(n = n_levels_fixeff[x], k = n/n_levels_fixeff[x]))
    if (length(variance_fixeff) != length(variables$covariates) & length(variance_fixeff) == 1){
      variances_fixed_effects <- rep(variance_fixeff, length(variables$covariates))
    } else {
      variances_fixed_effects <- variance_fixeff
    }
    values_fixed_effects <- lapply(seq_along(variables$covariates), function(x){
      pat_intercept <- c(0,rnorm(n_levels_fixeff[x]-1, 0, variances_fixed_effects[x])) # random intercept for each patient
      rep(pat_intercept, each = n/n_levels_fixeff[x])
    })
    covariates_list <- append(covariates_list,lapply(seq_along(variables$covariates), function(x){values_fixed_effects[[x]]}))
    for (i in seq_along(variables$covariates)) {
      out <- out %>% dplyr::mutate(!!variables$covariates[i] := values_fixed_effects[[i]])
    }
  }
  if (!is.null(variables$random_covariates) & (type == "glmer")){ # random effects covariates
    random_intercepts <- lapply(seq_along(n_levels_raneff), function(x) gl(n = n_levels_raneff[x], k = n/n_levels_raneff[x]))
    variances_random_intercepts <- rep(variance_raneff, length(n_levels_raneff))
    values_random_intercepts <- lapply(seq_along(random_intercepts), function(x){
      pat_intercept <- rnorm(n_levels_raneff[x], 0, variances_random_intercepts[x]) # random intercept for each patient
      rep(pat_intercept, each = n/n_levels_raneff[x])
    })
    list_random_covariates <- lapply(seq_along(variables$random_covariates), function(x){values_random_intercepts[[x]]})
    for (i in seq_along(variables$random_covariates)) {
      out <- out %>% dplyr::mutate(!!variables$random_covariates[i] := values_random_intercepts[[i]])
    }
  } else{
    list_random_covariates <- NULL
  }
  out_counts_ls <- purrr::map(seq_len(number_of_differential_clusters), function(i){
    linear_sum <- linear_sum(covariates_list = covariates_list,
                             betas = b[[i]],
                             random_covariates_list = list_random_covariates)
    if (type == "lm") {# linear regression
      Y <-  linear_sum + rnorm(n,0,error_variance) # True dependency
      Y_True <- linear_sum
    }else if (type %in% c("glm","glmer")) { # logistic regression
      Y <-  1/(1+exp(-linear_sum)) # True dependency
      Y_with_error <- Y + rnorm(length(Y),0,error_variance)
      Y_with_error <- ifelse(Y_with_error>1 | Y_with_error < 0, Y, Y_with_error)
      Y_size <-  unlist(purrr::map2(Y_with_error,size_tot, ~ rbinom(1, .y, .x )))
      Y <- Y_size/size_tot
      Y_True <- (1/(1+exp(-linear_sum)))
    }
    tmp_response_name <- ifelse(is.null(number_of_clusters), variables$response, paste0(variables$response,"_",i))
    y_tib <- tibble::tibble(!!tmp_response_name := as.vector(Y),
                            !!paste0(tmp_response_name,"_True") := as.vector(Y_True))
    if (!is.null(number_of_clusters)){
      stopifnot((number_of_clusters%%1==0) & (number_of_clusters > 2))
      if (type %in% c("glm","glmer")){
        remaining_y <- 1-Y # rest size
        Y_rand <- unlist(purrr::map(seq_along(remaining_y), function(i){ # random sizes for the remaining clusters without signal
          tmp <- runif(number_of_clusters-1,0,remaining_y[i])
          tmp_diff <- diff(sort(tmp))
          tmp_cumsum <- cumsum(tmp_diff)
          tmp_y_rand <- c(tmp_diff,remaining_y[i]-tmp_cumsum[number_of_clusters-2])
          return(tmp_y_rand)
        }))
        Y_rand <- matrix(Y_rand,nrow = number_of_clusters-1)
        Y_comp <- abs(rbind(Y,Y_rand))
      } else {
        Y_rand <- matrix(runif(n * (number_of_clusters-1), min(Y),max(Y)), ncol = n)
        Y_comp <- rbind(Y,Y_rand)
      }
      return(list(out = y_tib, d_counts = Y_comp))
    } else{
      return(list(out = y_tib))
    }

  })
  if (is.null(number_of_clusters)){
    return(dplyr::bind_cols(out,out_counts_ls[[1]]$out))
  } else {
    out <- dplyr::bind_cols(out,purrr::map(out_counts_ls, ~ .x$out))
    real_sizes <- unlist(purrr::map(out_counts_ls, function(counts){
      max(counts[[2]][1, ])
    }))
    wanted_sizes <- calculate_wanted_cluster_proportion(real_sizes)
    comb_counts <- purrr::map(seq_along(out_counts_ls), function(i){
      counts <- out_counts_ls[[i]][[2]]
      c <- (colSums(counts[-1, ]))/(wanted_sizes[i]-(1-colSums(counts[-1, ])))
      stopifnot(c>0)
      counts[-1, ] <- t(t(counts[-1, ])/c)
      return(counts)
    })
    Y_comp <- do.call(rbind,purrr::map(comb_counts, ~ .x))
    rownames(Y_comp) <- paste0("F_",seq_len(number_of_clusters*number_of_differential_clusters))
    rownames(Y_comp)[((0:(number_of_differential_clusters-1)))*number_of_clusters+1] <-
      gsub(pattern = "F",replacement = "T",
           x = rownames(Y_comp)[((0:(number_of_differential_clusters-1)))*number_of_clusters+1])
    colnames(Y_comp) <- seq_len(n)
    d_counts <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = round(t(out$size_tot*t(Y_comp)))),
      rowData = data.frame(cluster_id = rownames(Y_comp)))
out$size_tot <- colSums(SummarizedExperiment::assay(d_counts))
    return(list(out = out, d_counts = d_counts))
  }
}
