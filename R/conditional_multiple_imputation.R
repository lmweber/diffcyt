#' Conditional multiple imputation
#'
#' Imputes censored values according to \href{https://www.researchgate.net/publication/319246304_Improved_conditional_imputation_for_linear_regression_with_a_randomly_censored_predictor}{Atem et al. 2017}
#'
#' @param data 'data.frame'
#' @param formula the formula for fitting the regression model with a special
#'  syntax for the censored covariate : e.g. 'y~Surv(x,I)' means 'y~x' with 'x' being
#'  censored and 'I' the censoring indicator (0=censored,1=observed).
#' @param regression_type function. The regression type to be used, lm for linear
#' regression, glm for general linear regression, glmer for generalized
#' linear mixed-effects models. Default: lm
#' @param repetitions number of repetitions for multiple imputation. Default: 10
#' @param method_est which method should be used in the imputation step. One of
#'  'cc', 'pmm', 'mrl', 'rs', 'km'. See details.
#' @param weights Weights to be used in fitting the regression model. Default = NULL
#' @param contrasts Contrast vector to be used in testing the regression model. Default = NULL
#' @param family The family to be used in the regression model. Default = "binomial". Omited if linear model is used.
#'  for the censored values.
#' @param id name of column containing id of sample
#' @param verbose Logical.
#' @param n_obs_min minimum number of observed events needed. default = 2.
#'  if lower than this value will throw an error.
#' @param  ... additional arguments passed to testing
#'
#' @details Possible methods in 'methods_est' are:
#' \describe{
#'   \item{'cc'}{complete case}
#'   \item{'pmm'}{predictive mean matching, treating censored values as missing}
#'   \item{'mrl'}{Mean Residual Life (Conditional single imputation from \href{https://www.researchgate.net/publication/319246304_Improved_conditional_imputation_for_linear_regression_with_a_randomly_censored_predictor}{Atem et al. 2017})}
#'   \item{'rs'}{Risk Set imputation}
#'   \item{'km'}{Kaplan Meier imputation}
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom lme4 glmer
#' @export
#' @examples
#'  lm_formula <- formula(Y ~ Surv(X,I) + Z)
#'  data <- simulate_data(100, lm_formula, type = "lm", n_levels_fixeff=2)
#'  cmi_out <- conditional_multiple_imputation(data,lm_formula)
conditional_multiple_imputation <-
  function(data,
           formula,
           regression_type = c("lm", "glm", "glmer"),
           repetitions = 10,
           method_est = c("mrl", "rs", "km", "cc","pmm"),
           weights = NULL,
           contrasts = NULL,
           family = "binomial",
           id = NULL,
           verbose = FALSE,
           n_obs_min = 2,
           ...) {
  data <- dplyr::as_tibble(data)
  n <- dim(data)[1]
  variables_formula <- extract_variables_from_formula(formula)
  formula <- create_glmm_formula(formula)
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
  method_est <- match.arg(method_est)
  regression_type <- match.arg(regression_type)
  if (!is.null(weights)) {
    data$weights <- weights
  }
  data <-
    data_processing_for_imputation(data,
                             censored_variable,
                             censoring_indicator,
                             response,
                             covariates ,
                             id)
  if (sum(data[[censoring_indicator]]) <= n_obs_min){
    stop(paste0("Not enough observed events, consider decreasing 'n_obs_min (currently = ",n_obs_min," )"),call. = FALSE)
  }
  # complete case
  if (method_est == "cc"){
    return(
      complete_case(
        data = data,
        censored_variable = censored_variable,
        censoring_indicator = censoring_indicator,
        formula = formula,
        regression_type = regression_type,
        weights = "weights",
        family = family
      )
    )
    # predictive mean matching (treating censored values as missing)
  } else if (method_est == "pmm"){
    return(
      predictive_mean_matching(
        data = data,
        censored_variable = censored_variable,
        censoring_indicator = censoring_indicator,
        variables_for_imputation = c(response,covariates,variables_formula[["random_covariates"]]),
        formula = formula,
        regression_type = regression_type,
        repetitions = repetitions,
        weights = "weights",
        family = family
      )
    )
  }
  formula <- sub(censored_variable,est_name,Reduce(paste, deparse(formula)))

  # iterate
  all_csi_out <- lapply(seq_len(repetitions), function(i){
      # if CSI, no random sampling
      if (repetitions == 1 | method_est %in% c("rs","km")) {
        randsam <- seq_len(n)
      } else {
        # random sample same length
        randsam <- sample(seq_len(n), size = n, replace = TRUE)
      }
      return(
        conditional_single_imputation(
          data = data[randsam, ],
          censored_variable = censored_variable,
          censoring_indicator = censoring_indicator,
          response = response,
          covariates = covariates,
          id = id,
          method_est = method_est,
          verbose = verbose
        )
      )
    })

  data_cmi <-
    conditional_multiple_imputation_testing(
      data = data,
      imputed_datasets = all_csi_out,
      censored_variable = censored_variable,
      censoring_indicator = censoring_indicator,
      response = response,
      covariates = covariates,
      id = id,
      formula = formula,
      regression_type = regression_type,
      repetitions = repetitions,
      method_est = method_est,
      verbose = verbose,
      weights = weights,
      contrasts = contrasts,
      family = family,
      ...
    )
  return(data_cmi)
}

#' perform fitting of imputed datasets
#' @param data 'data.frame'
#' @param imputed_datasets all the imputed datasets for which testing is to be
#'  performed.
#' @param censored_variable name of column containing censored data
#' @param censoring_indicator name of column containing indication if observed
#'   (1) or censored (0) value in column 'censored_variable'
#' @param response name of column containing the response (can be 'NULL')
#' @param covariates name(s) of column(s) containing the covariate that influences
#'  censoring
#' @inheritParams conditional_multiple_imputation
#' @importFrom magrittr %>%
#' @importFrom lme4 glmer
#' @importFrom  stats na.omit
conditional_multiple_imputation_testing <-
  function(data,
           imputed_datasets,
           censored_variable,
           censoring_indicator,
           response = NULL,
           covariates = NULL,
           id = NULL,
           formula = NULL,
           regression_type = c("lm", "glm", "glmer"),
           repetitions = 10,
           method_est = c("mrl", "rs", "km"),
           verbose = FALSE,
           weights = NULL,
           contrasts = NULL,
           family = "binomial",
           ...) {

    censored_variable <- censored_variable
    expr <- rlang::enquo(censored_variable)
    est_name <- paste0(rlang::quo_name(expr),"_est")
    if (is.null(id)) {
      id_was_null <- TRUE
      data$id <- seq_along(data[[censored_variable]])
      id <- "id"
    } else{
      id_was_null <- FALSE
    }
    regression_type <- match.arg(regression_type)

    fits <- lapply(imputed_datasets, function(csi_out) {
      if (id_was_null) {
        csi_out$id <- data$id
      }
      max_est_name <- max(na.omit(csi_out[[est_name]]))
      cond_2 <- ifelse(is.na(max_est_name),
                       FALSE,
                       abs(max_est_name) < 1000 * abs(max(csi_out[[censored_variable]])))
      if (!any(is.na(csi_out[[censored_variable]])) &
          cond_2 &
          (sum(csi_out[[censoring_indicator]]) > 2)) {
        return(tryCatch({
          args <- args_for_fitting(csi_out, formula, regression_type, family)
          do.call(regression_type, args)
        }, error = function(e) {
          csi_out_red <- dplyr::filter(csi_out, is.finite(csi_out[[est_name]]))
          if (verbose)
            message(paste("NaN, NA of Inf present, removing",
                dim(csi_out)[1] - dim(csi_out_red)[1], "/", dim(csi_out)[1],
                "values for fitting"))
          args <- args_for_fitting(csi_out_red, formula, regression_type, family)
          do.call(regression_type, args)
        }))
      }
    })

    # prepare output matrix
    ls_EstX <- dplyr::tibble(!!id := data[[id]])
    # fits <- list()

    for (i in seq_len(repetitions)) {
      if (id_was_null) {
        imputed_datasets[[i]]$id <- data$id
      }
      csi_out <- imputed_datasets[[i]]
    ls_EstX <- suppressMessages(
      csi_out %>%
        dplyr::distinct() %>%
        dplyr::select(!!est_name, !!id) %>%
        dplyr::group_by(!!dplyr::sym(id)) %>%
        dplyr::summarise_at(est_name, list(mean = mean)) %>%
        dplyr::right_join(ls_EstX) %>%
        dplyr::rename(!!paste0(est_name, "_", i) := mean)
    )
    }
    fits <- fits[unlist(purrr::map(fits, ~ !is.null(.x)))]
    betas <- get_betas_from_multiple_fits(fits,regression_type)
    betasVar <- get_variance_of_betas_from_multiple_fits(fits)
    ls_EstX <- ls_EstX %>%
      dplyr::arrange(!!dplyr::sym(id)) %>%
      dplyr::select(!!id, dplyr::starts_with(est_name))

    # mean imputed value
    EstXMean <-
      apply(ls_EstX %>% dplyr::select(-!!id), 1, function(x)
        mean(na.omit(x)))
    data <- dplyr::arrange(data, !!dplyr::sym(id))
    data <-
      tryCatch(
        data %>% dplyr::mutate(!!est_name := EstXMean),
        error = function(e) data)
    data <- data[ , !colnames(data)=="rank"]
    data_cmi <- list(
      data = data,
      betasMean = colMeans(betas),
      betasVar = colMeans(betasVar),
      metadata = list(
        repetitions = repetitions,
        betas = betas,
        vars = betasVar#,
        # ls_EstX = ls_EstX
      ),
      fits = fits
    )
    return(data_cmi = data_cmi)
  }
