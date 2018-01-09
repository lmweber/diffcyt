#' Create model formula and data frame of variables
#' 
#' Create model formula and corresponding data frame of variables for model fitting
#' 
#' Creates a model formula and corresponding data frame of variables specifying the models
#' to be fitted. (Alternatively, \code{\link{createDesignMatrix}} can be used to generate
#' a design matrix instead of a model formula.)
#' 
#' The output is a list containing the model formula and corresponding data frame of
#' variables (one column per formula term). These can then be provided to differential
#' testing functions that require a model formula, together with the main data object and
#' contrast matrix.
#' 
#' Depending on the experimental design, there may also be block IDs (e.g. patient IDs in
#' a paired design), batch effects, and continuous covariates. These can be included in
#' the model formula by providing the optional arguments \code{block_IDs},
#' \code{batch_IDs}, and \code{covariates}.
#' 
#' Block IDs can be included as either fixed effect terms or random intercept terms. This
#' choice is controlled with the \code{block_IDs_type} argument. (See arguments for more
#' details.)
#' 
#' In addition, sample IDs can be included as random intercept terms by providing the
#' \code{sample_IDs} argument, in order to account for the overdispersion typically seen
#' in high-dimensional cytometry data. This is known as an 'observation-level random
#' effect'; see Nowicka et al. (2017), \emph{F1000Research}, for more details.
#' 
#' 
#' @param group_IDs Vector or factor of group membership labels for each sample (e.g.
#'   diseased vs. healthy, or treated vs. untreated). These are the parameters of interest
#'   in the model. Vectors are converted to factors internally. The first level of the
#'   factor (e.g. healthy) or first entry of the vector will be used as the reference
#'   level for differential testing. To re-order factor levels, use
#'   \code{\link[stats]{relevel}} or \code{\link[base]{factor}}.
#' 
#' @param block_IDs (Optional) Vector or factor of block IDs (e.g. patient IDs for paired
#'   designs, such as one diseased and one healthy sample per patient). Block IDs can be
#'   included as either fixed effect terms or random intercept terms; this choice is
#'   specified using the argument \code{block_IDs_type}. Note that some testing methods
#'   (\code{\link{testDA_limma}} and \code{\link{testDS_med}}) require block IDs to be
#'   provided directly to the testing function, if they are to be included as random
#'   intercepts.
#' 
#' @param block_IDs_type (Optional) Whether block IDs should be included in the model
#'   formula as fixed effect terms or random intercept terms. Options are 'fixed' and
#'   'random'. Default = 'fixed'. Note that some testing methods
#'   (\code{\link{testDA_limma}} and \code{\link{testDS_med}}) require block IDs to be
#'   provided directly to the testing function, if they are to be included as random
#'   intercepts.
#' 
#' @param batch_IDs (Optional) Vector or factor of batch effect IDs. Batch effect IDs are
#'   included as fixed effect terms in the model formula. This allows batch effects to be
#'   estimated during model fitting, and taken into account during inference on the group
#'   ID parameters of interest.
#' 
#' @param covariates (Optional) Numeric matrix of continuous covariates (one column per
#'   covariate). Covariates are included as fixed effect terms in the model formula. This
#'   allows their effects to be estimated during model fitting, and taken into account
#'   during inference on the group ID parameters of interest.
#' 
#' @param sample_IDs (Optional) Vector or factor of sample IDs. Sample IDs are included as
#'   random intercept terms in the model formula, to account for the overdispersion
#'   typically seen in high-dimensional cytometry data. This is known as an
#'   'observation-level random effect'; for more details see Nowicka et al. (2017),
#'   \emph{F1000Research}.
#' 
#' 
#' @return Returns a list with two elements:
#' \itemize{
#' \item \code{formula}: model formula
#' \item \code{data}: data frame of variables corresponding to the model formula
#' }
#' 
#' 
#' @importFrom stats as.formula
#' @importFrom methods is
#' 
#' @export
#' 
#' @examples
#' # A full workflow example demonstrating the use of each function in the 'diffcyt'
#' # pipeline on an experimental data set is available in the package vignette.
#' 
createFormula <- function(group_IDs, 
                          block_IDs = NULL, block_IDs_type = c("fixed", "random"), 
                          batch_IDs = NULL, covariates = NULL, 
                          sample_IDs = NULL) {
  
  if (!is.factor(group_IDs)) {
    group_IDs <- factor(group_IDs, levels = unique(group_IDs))
  }
  if (!is.null(block_IDs) & !is.factor(block_IDs)) {
    block_IDs <- factor(block_IDs, levels = unique(block_IDs))
  }
  if (!is.null(batch_IDs) & !is.factor(batch_IDs)) {
    batch_IDs <- factor(batch_IDs, levels = unique(batch_IDs))
  }
  if (!is.null(covariates) & !(is.matrix(covariates) & is.numeric(covariates))) {
    stop("'covariates' must be provided as a numeric matrix, with one column for each covariate")
  }
  if (!is.null(sample_IDs) & !is.factor(sample_IDs)) {
    sample_IDs <- factor(sample_IDs, levels = unique(sample_IDs))
  }
  
  # create formula
  
  # fixed effects
  terms <- c("group_IDs", "block_IDs", "batch_IDs", "covariates")
  nulls <- c(is.null(group_IDs), 
             is.null(block_IDs) | block_IDs_type == "random", 
             is.null(batch_IDs), 
             is.null(covariates))
  
  formula_chr <- paste("y ~", paste(terms[!nulls], collapse = " + "))
  
  # random effects
  if (!is.null(block_IDs) & block_IDs_type == "random") {
    formula_chr <- paste0(formula_chr, " + (1 | block_IDs)")
  }
  if (!is.null(sample_IDs)) {
    formula_chr <- paste0(formula_chr, " + (1 | sample_IDs)")
  }
  
  formula <- as.formula(formula_chr)
  
  # create data frame of variables (in correct order)
  data <- data.frame(group_IDs)
  
  if (!is.null(block_IDs) & block_IDs_type == "fixed") data <- cbind(data, block_IDs)
  if (!is.null(batch_IDs)) data <- cbind(data, batch_IDs)
  if (!is.null(covariates)) data <- cbind(data, covariates)
  if (!is.null(block_IDs) & block_IDs_type == "random") data <- cbind(data, block_IDs)
  if (!is.null(sample_IDs)) data <- cbind(data, sample_IDs)
  
  # return as list
  list(formula = formula, data = data)
}


