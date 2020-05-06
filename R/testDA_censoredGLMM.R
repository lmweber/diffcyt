#' Test for differential abundance: method 'censcyt-DA-GLMM'
#'
#' @param d_counts \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'  object containing cluster cell counts, from \code{\link[diffcyt]{calcCounts}}.
#' @param formula See \code{\link[diffcyt]{createFormula}}
#'  for details but make sure to specify the 'formula$formula' as follows: The
#'  censored predictor is encoded as Surv(censored_variable, censoring_indicator)'.
#'  e.g. y~Surv(x,I)
#' @param contrast Contrast matrix, created with
#'  \code{\link[diffcyt]{createContrast}}.
#' @param m number of repetition for multiple imputation. default is m=10.
#' @param method_est which method should be used in the imputation step. One of
#'  'cc', 'pmm', 'mrl', 'rs', 'km'. See details.
#' @param min_cells positive Integer. The minimum number of cells a population
#'  needs to have to be included. Default = 3.
#' @param min_samples postive Integer. The minimum number of samples to still
#'  fit a model. Default = 3.
#' @param normalize Whether to include optional normalization factors to adjust for
#'   composition effects (see details). Default = FALSE.
#' 
#' @param norm_factors Normalization factors to use, if \code{normalize = TRUE}. Default =
#'   \code{"TMM"}, in which case normalization factors are calculated automatically using
#'   the 'trimmed mean of M-values' (TMM) method from the \code{edgeR} package.
#'   Alternatively, a vector of values can be provided (the values should multiply to 1).
#' 
#' @param BPPARAM specify parallelization option as one of 
#'  \code{\link[BiocParallel]{BiocParallelParam}}.
#'  e.g. \code{\link[BiocParallel]{MulticoreParam}}(workers=2) for parallelization 
#'  with two cores. Default = \code{\link[BiocParallel]{SerialParam}} 
#'  (no parallelization). Parallelization works only if 'BiocParallel' is available.
#' @param verbose Logical.
#' 
#' @details Possible methods in 'methods_est' are:
#' \describe{
#'   \item{'cc'}{complete case, removing incomlete samples}
#'   \item{'pmm'}{predictive mean matching, treating censored values as missing}
#'   \item{'mrl'}{Mean Residual Life (Conditional single imputation from \href{https://www.researchgate.net/publication/319246304_Improved_conditional_imputation_for_linear_regression_with_a_randomly_censored_predictor}{Atem et al. 2017})}
#'   \item{'rs'}{Risk Set imputation}
#'   \item{'km'}{Kaplan Meier imputation}
#' }
#'
#' @export
#' @examples 
#' tmp_formula <- formula(y~Surv(X,I)+z+(1|r))
#' data_sim <- simulate_data(
#'  n = 20,
#'  formula = tmp_formula,
#'  n_levels_fixeff = 2,
#'  type = "glmer",
#'  b = list(b=c(-5,-2,0.2)),
#'  number_of_clusters = 10,
#'  number_of_differential_clusters = 1)
#'
#' d_counts <- data_sim[["d_counts"]]
#' data_sim <- data_sim[["out"]]
#' da_formula <- list(formula = tmp_formula,
#'                   data = dplyr::select(data_sim,"X","I","z","r"),
#'                   random_terms = TRUE)
#' contrast <- diffcyt::createContrast(c(0, 1, 0))
#'
#' outs <- testDA_censoredGLMM(d_counts = d_counts, formula = da_formula,
#'                            contrast = contrast, method_est = "km",
#'                            verbose = FALSE, m = 10)
#'
testDA_censoredGLMM <- function(d_counts, formula, contrast, m = 10,
                                method_est = c("km","rs","mrl","cc","pmm"),
                                min_cells = 3,
                                min_samples = 3,
                                normalize = FALSE, 
                                norm_factors = "TMM",
                                BPPARAM=BiocParallel::SerialParam(),
                                verbose = FALSE
                                )
{
  method_est <- match.arg(method_est)
  BPPARAM <- if(requireNamespace("BiocParallel")){BPPARAM} else{NULL}
  # variable names from the given formula
  cmi_input <- extract_variables_from_formula(formula$formula)
  # create formula for fitting
  formula_glmm <- create_glmm_formula(formula$formula)
  
  counts <- SummarizedExperiment::assays(d_counts)[["counts"]]
  cluster_id <- SummarizedExperiment::rowData(d_counts)$cluster_id
  # only keep counts with more than the minimum number of cells
  counts_to_keep <- counts >= min_cells
  # only keep clusters with more than the minimum number of samples
  rows_to_keep <- apply(counts_to_keep, 1, function(r) sum(r) >= min_samples)
  # subset counts and cluster_id's
  counts <- counts[rows_to_keep, , drop = FALSE]
  cluster_id <- cluster_id[rows_to_keep]

    # normalization factors
  if (normalize & norm_factors == "TMM") {
    norm_factors <- calcNormFactors(counts, method = "TMM")
  }
  
  if (normalize) {
    weights <- colSums(counts) / norm_factors
  } else {
    weights <- colSums(counts)
  }
  
  if (ncol(contrast) == 1 & nrow(contrast) > 1) {
    contrast <- t(contrast)
  }

  if (verbose) message(paste(sum(formula$data[[cmi_input[["censoring_indicator"]]]] == 0),
                           "of", dim(formula$data)[1], "values are censored"))
  # start fitting by iterating through the clusters, use normal lapply if 
  # 'BiocParallel' isn't install otherwise use bplapply
  p_vals_ls <- maybe_parallel_lapply(seq_along(cluster_id), 
                                     BPPARAM=BPPARAM, 
                                     function(i) {
    # data preparation
    y <- counts[i, ]/weights
    data_i <- cbind(y, weights, formula$data)
    colnames(data_i)[c(1,2)] <- c(cmi_input$response,"weights")
    
    # do conditional multiple imputation
    if (method_est %in% c("mrl","rs","km","pmm")){
      out_test <- tryCatch(suppressMessages(suppressWarnings(
        conditional_multiple_imputation(
          data = data_i,
          formula = formula$formula,
          repetitions = m,
          method_est = method_est,
          regression_type = "glmer",
          family = "binomial",
          verbose = verbose,
          weights = weights,
          contrasts = contrast
        ))),
      error=function(e) NA)
      # pooling of results from multiple imputation
      p_val <- tryCatch(summary(mice::pool(out_test$fits))$p.value[ which(contrast == 1)],
                        error=function(e) NA)
      if (verbose) {
        hmi <-  tryCatch(how_many_imps(out_test), error = function(x) {NA})
        message(paste("min suggested number of imputations ('m')", ceiling(hmi)))
      }
    } # complete case fitting and testing
    else if (method_est == "cc"){
      p_val <- tryCatch({
      fit <- suppressMessages(suppressWarnings(complete_case(
        data = data_i,censored_variable = cmi_input[["censored_variable"]],
        censoring_indicator = cmi_input[["censoring_indicator"]],
        formula = formula_glmm,regression_type = "glmer",
        weights = "weights",family = "binomial")$fits))
      test <- multcomp::glht(fit, contrast)
      summary(test)$test$pvalues
      },error=function(e) NA)
    }
    return(p_val)
  })
  p_vals <- unlist(p_vals_ls)
  # fdr correction
  p_adj <- p.adjust(p_vals, method = "fdr")
  stopifnot(length(p_vals) == length(p_adj))
  
  row_data <- data.frame(cluster_id = as.character(cluster_id),
                         p_val = p_vals,
                         p_adj = p_adj,
                         stringsAsFactors = FALSE)
  # return NA if cluster has been excluded from testing because of too few observations
  row_data <- suppressMessages(
    dplyr::right_join(row_data,
                      data.frame(cluster_id = as.character(SummarizedExperiment::rowData(d_counts)$cluster_id),stringsAsFactors = FALSE)))
  
  res <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = SummarizedExperiment::assay(d_counts)),
    rowData = row_data)
  
  # return normalization factors in 'metadata'
  if (normalize) {
    metadata(res)$norm_factors <- norm_factors
  }
  return(res)
}

# if namespace 'BiocParallel' is available run 'BiocParallel::bplapply', otherwise
# run lapply
maybe_parallel_lapply <- function(X, FUN, BPPARAM) {
  if (requireNamespace("BiocParallel", quietly = TRUE)) {
    BiocParallel::bplapply(X, FUN, BPPARAM=BPPARAM)
  } else {
    lapply(X, FUN)
  }
}