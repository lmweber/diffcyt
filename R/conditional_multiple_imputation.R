#' Conditional multiple imputation
#'
#' First two steps for multiple imputation for censored covariates. Returns
#' regression fits in a list that can be combined using \code{\link[mice]{pool}}().
#' 
#' 
#'
#' @param data 'data.frame'
#' 
#' @param formula the formula for fitting the regression model with a special
#'  syntax for the censored covariate : e.g. 'y~Surv(x,I)' means 'y~x' with 'x' being
#'  censored and 'I' the event indicator (0=censored,1=observed).
#'  
#' @param regression_type function. The regression type to be used, lm for linear
#' regression, glm for general linear regression, glmer for generalized
#' linear mixed-effects models. Default: lm
#' 
#' @param mi_reps number of repetitions for multiple imputation. Default: 10
#' 
#' @param imputation_method which method should be used in the imputation step. One of
#'  'km','km_exp','km_wei','km_os', 'rs', 'mrl', 'cc', 'pmm'. See details. default = 'km'.
#'  
#' @param weights Weights to be used in fitting the regression model. Default = NULL
#' 
#' @param contrasts Contrast vector to be used in testing the regression model. 
#' Default = NULL
#' 
#' @param family The family to be used in the regression model. Default = "binomial". 
#' Omitted if linear model is used.
#' 
#' @param id name of column containing id of sample
#' 
#' @param verbose Logical.
#' 
#' @param n_obs_min minimum number of observed events needed. default = 2.
#'  if lower than this value will throw an error.
#' 
#' 
#' @references {
#'  A Comparison of Several Methods of Estimating the Survival Function When 
#'  There is Extreme Right Censoring (M. L. Moeschberger and John P. Klein, 1985)
#'  }
#' @details Possible methods in 'imputation_method' are:
#' \describe{
#'   \item{'km'}{Kaplan Meier imputation is similar to 'rs' (Risk set imputation) 
#'               but the random draw is according to the survival function of
#'               the respective risk set.}
#'   \item{'km_exp'}{The same as 'km' but if the largest value is censored the 
#'              tail of the survival function is modeled as an exponential 
#'              distribution where the rate parameter is obtained by fixing
#'              the distribution to the last observed value. 
#'              See (Moeschberger and Klein, 1985).}
#'   \item{'km_wei'}{The same as 'km' but if the largest value is censored the 
#'              tail of the survival function is modeled as an weibull 
#'              distribution where the parameters are obtained by MLE fitting on
#'              the whole data. See (Moeschberger and Klein, 1985).}
#'   \item{'km_os'}{The same as 'km' but if the largest value is censored the 
#'              tail of the survival function is modeled by order statistics. 
#'              See (Moeschberger and Klein, 1985).}
#'   \item{'rs'}{Risk Set imputation replaces the censored values with a random
#'               draw from the risk set of the respective censored value.}
#'   \item{'mrl'}{Mean Residual Life (Conditional single imputation from 
#'                \href{https://www.researchgate.net/publication/319246304_Improved_conditional_imputation_for_linear_regression_with_a_randomly_censored_predictor}{Atem et al. 2017})
#'                is a multiple imputation procedure that bootstraps the data and
#'                imputes the censored values by replacing them with their 
#'                respective mean residual life.}
#'   \item{'cc'}{complete case (listwise deletion) analysis removes incomlete samples.}
#'   \item{'pmm'}{predictive mean matching treats censored values as missing and
#'                uses predictive mean matching method from \code{\link[mice]{mice}}.}
#' }
#'  
#'  
#' @return A list with five elements:
#' \describe{
#'   \item{'data'}{The input data frame}
#'   \item{'betasMean'}{the mean regression coefficients}
#'   \item{'betasVar'}{the variances of the mean regression coefficients}
#'   \item{'metadata'}{a list of three elements: \describe{
#'     \item{'mi_reps'}{number of repetitions in multiple imputation}
#'     \item{'betas'}{all regression coefficients}
#'     \item{'vars'}{the variances of the regression coefficients}
#'   }}
#'   \item{'fits'}{list with all regression fits}
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom lme4 glmer
#' @importFrom rlang :=
#' @export
#' @examples
#'  lm_formula <- formula(Y ~ Surv(X,I) + Z)
#'  data <- simulate_singlecluster(100, lm_formula, type = "lm", n_levels_fixeff=2)
#'  cmi_out <- conditional_multiple_imputation(data,lm_formula)
#'  comb_out <- mice::pool(cmi_out$fits)
#'  pvals <- summary(comb_out)$p.value
conditional_multiple_imputation <-
  function(data,
           formula,
           regression_type = c("lm", "glm", "glmer"),
           mi_reps = 10,
           imputation_method = c("km","km_exp","km_wei","km_os", "rs", "mrl", "cc", "pmm"),
           weights = NULL,
           contrasts = NULL,
           family = "binomial",
           id = NULL,
           verbose = FALSE,
           n_obs_min = 2
           ) {
  data <- dplyr::as_tibble(data)
  n <- dim(data)[1]
  variables_formula <- extract_variables_from_formula(formula)
  formula_uncens <- create_glmm_formula(formula)
  censored_variable <- variables_formula[["censored_variable"]]
  censoring_indicator <- variables_formula[["censoring_indicator"]]
  response <- variables_formula[["response"]]
  covariates <- variables_formula[["covariates"]]

  
  expr <- rlang::enquo(censored_variable)
  est_name <- paste0(rlang::quo_name(expr),"_est")
  if (is.null(id)) {
    data$id <- seq_along(data[[censored_variable]])
    id <- "id"
  }
  if (is.numeric(id) & !is.null(id)) {
    id <- colnames(data)[[id]]
  }
  binarise_covariate <- stringr::str_detect(imputation_method,"_bin")
  imputation_method <- ifelse(binarise_covariate,stringr::str_split(imputation_method,"_bin")[[1]][1],imputation_method)
  imputation_method <- match.arg(imputation_method)
  regression_type <- match.arg(regression_type)
  
  if (is.character(weights)) {
    stopifnot(all(weights %in% colnames(data)))
  } else if (is.numeric(weights) | is.logical(weights)) {
    weights <- colnames(data)[weights]
  }


  # check that data is in correct format, some type conversions, but keep covariates = NULL to not convert to factors
  data <- data_processing_for_imputation(data, 
                                         censored_variable, 
                                         censoring_indicator, 
                                         covariates = NULL , 
                                         id)
  
  if (sum(data[[censoring_indicator]]) <= n_obs_min){
    stop(paste0("Not enough observed events, consider decreasing 'n_obs_min (currently = ",n_obs_min," )"),call. = FALSE)
  }
  # complete case
  if (imputation_method == "cc"){
    return(
      complete_case(
        data = data,
        censored_variable = censored_variable,
        censoring_indicator = censoring_indicator,
        formula = formula_uncens,
        regression_type = regression_type,
        weights = weights,
        family = family,
        binarise_covariate = binarise_covariate
      )
    )
    # predictive mean matching (treating censored values as missing)
  } else if (imputation_method == "pmm"){
    return(
      predictive_mean_matching(
        data = data,
        censored_variable = censored_variable,
        censoring_indicator = censoring_indicator,
        variables_for_imputation = c(response,covariates,variables_formula[["random_covariates"]]),
        formula = formula_uncens,
        regression_type = regression_type,
        mi_reps = mi_reps,
        weights = weights,
        family = family
      )
    )
  }

  data_spread <- create_two_level_factor_data(data, covariates)
  covariates_spread <- create_two_level_factor_covariates(data, formula)
  # iterate
  if (mi_reps > 1 & imputation_method %in% c("km","km_exp","km_wei","km_os","rs")){
    all_csi_out <- estimate_cens_vars(data_spread,censored_variable,censoring_indicator,covariates,id,imputation_method,mi_reps)
  }
  else{
  all_csi_out <- lapply(seq_len(mi_reps), function(i){
      # if CSI, no random sampling
      if (mi_reps == 1 | imputation_method %in% c("rs","km")) {
        randsam <- seq_len(n)
      } else {
        # random sample same length
        randsam <- sample(seq_len(n), size = n, replace = TRUE)
      }

      return(
        conditional_single_imputation(
          data = data_spread[randsam, ],
          censored_variable = censored_variable,
          censoring_indicator = censoring_indicator,
          response = response,
          covariates = covariates_spread,
          id = id,
          imputation_method = imputation_method,
          verbose = verbose
        )
      )
    })
  }
  formula_uncens_est_name <- sub(censored_variable,est_name,Reduce(paste, deparse(formula_uncens)))
  
  data_cmi <-
    conditional_multiple_imputation_fitting(
      data = data,
      imputed_datasets = all_csi_out,
      censored_variable = censored_variable,
      censoring_indicator = censoring_indicator,
      response = response,
      covariates = covariates,
      id = id,
      formula = formula_uncens_est_name,
      regression_type = regression_type,
      mi_reps = mi_reps,
      verbose = verbose,
      weights = weights,
      contrasts = contrasts,
      family = family,
      binarise_covariate = binarise_covariate
    )
  return(data_cmi)
}

# perform fitting of imputed datasets
# @param data 'data.frame'
# @param imputed_datasets all the imputed datasets for which testing is to be
#  performed.
# @param censored_variable name of column containing censored data
# @param censoring_indicator name of column containing indication if observed
#   (1) or censored (0) value in column 'censored_variable'
# @param response name of column containing the response (can be 'NULL')
# @param covariates name(s) of column(s) containing the covariate that influences
#  censoring
# @inheritParams conditional_multiple_imputation
#' @importFrom magrittr %>%
#' @importFrom lme4 glmer
#' @importFrom  stats na.omit
#' @importFrom rlang :=
conditional_multiple_imputation_fitting <-
  function(data,
           imputed_datasets,
           censored_variable,
           censoring_indicator,
           response = NULL,
           covariates = NULL,
           id = NULL,
           formula = NULL,
           regression_type = c("lm", "glm", "glmer"),
           mi_reps = 10,
           verbose = FALSE,
           weights = NULL,
           contrasts = NULL,
           family = "binomial",
           binarise_covariate = FALSE
           ) {
    # name of imputed variable
    est_name <- paste0(censored_variable,"_est")
    regression_type <- match.arg(regression_type)
    if (is.character(weights)) {
      stopifnot(all(weights %in% colnames(data)))
    } else if (is.numeric(weights) | is.logical(weights)) {
      weights <- colnames(data)[weights]
    }
    censored_bool <- data[[censoring_indicator]] == 0
    # do fitting of all imputed data sets
    fits <- lapply(imputed_datasets, function(csi_out) {
      
      if(is(csi_out,"vector")){
        if(sum(censored_bool) != length(csi_out)){
          censored_bool[length(censored_bool)] <- FALSE
          stopifnot(sum(censored_bool)==length(csi_out))
        }
        est_col <- data[[censored_variable]]
        est_col[censored_bool] <- csi_out
        csi_out <- cbind(data,est_col)
        colnames(csi_out)[dim(csi_out)[2]] <- est_name
      }
      # some checks of imputed values, only do fitting if imputation gave reasonable results
      max_est_name <- max(na.omit(csi_out[[est_name]]))
      cond_2 <- ifelse(is.na(max_est_name),
                       FALSE,
                       abs(max_est_name) < 1000 * abs(max(csi_out[[censored_variable]])))
      if (!any(is.na(csi_out[[censored_variable]])) &
          cond_2 &
          (sum(csi_out[[censoring_indicator]]) > 2)) {
        if (binarise_covariate){
          csi_out[[est_name]] <- binarised_covariate(csi_out[[est_name]])
        }
        # fitting
        return(tryCatch({
          args <- args_for_fitting(csi_out, formula, regression_type, family, weights)
          do.call(regression_type, args)
          
        # remove not finite imputed values and try fitting again if error
        }, error = function(e) {
          csi_out_red <- dplyr::filter(csi_out, is.finite(csi_out[[est_name]]))
          if (verbose)
            message(paste("NaN, NA of Inf present, removing",
                dim(csi_out)[1] - dim(csi_out_red)[1], "/", dim(csi_out)[1],
                "values for fitting"))
          args <- args_for_fitting(csi_out_red, formula, regression_type, family, weights)
          do.call(regression_type, args)
        }))
      }
    })
    # remove empty fits
    fits <- fits[unlist(purrr::map(fits, ~ !is.null(.x)))]
    # get regression coefficients
    betas <- get_betas_from_multiple_fits(fits,regression_type)
    # get variances of regression coefficients
    betasVar <- get_variance_of_betas_from_multiple_fits(fits)

    data_cmi <- list(
      data = data,
      betasMean = colMeans(betas),
      betasVar = colMeans(betasVar),
      metadata = list(
        mi_reps = mi_reps,
        betas = betas,
        vars = betasVar
      ),
      fits = fits
    )
    return(data_cmi = data_cmi)
  }
