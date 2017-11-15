#' Create model formula and data frame of variables
#' 
#' Create model formula and corresponding data frame of variables for model fitting
#' 
#' Creates a model formula and corresponding data frame of variables specifying the models
#' to be fitted. Similar to \code{\link{createDesignMatrix}}, but returns a model formula
#' instead of a design matrix.
#' 
#' The output is a list containing the model formula and corresponding data frame (one
#' column per formula term). These can then be provided to differential testing functions
#' that require a model formula, together with the main data object and contrast matrix.
#' 
#' 
#' @param group_IDs Vector or factor of group membership labels for each sample (e.g.
#'   diseased vs. healthy, or treated vs. untreated). These are the parameters of interest
#'   in the model. Vectors are converted to factors internally. The first level of the
#'   factor (e.g. healthy) or first entry of the vector will be used as the reference
#'   level for differential testing. To re-order factor levels, use
#'   \code{\link[stats]{relevel}} or \code{\link[base]{factor}}.
#' 
#' @param batch_IDs (Optional) Vector or factor of batch IDs. Batch IDs are included as
#'   fixed effect terms in the model formula. This allows batch effects to be estimated
#'   during model fitting and taken into account during inference on the parameters of
#'   interest.
#' 
#' @param covariates (Optional) Numeric matrix of continuous covariates. Covariates are
#'   included as fixed effect terms in the model formula. This allows their effects to be
#'   estimated during model fitting and taken into account during inference on the
#'   parameters of interest.
#' 
#' @param block_IDs (Optional) Vector or factor of block IDs, e.g. for paired designs
#'   (e.g. one diseased and one healthy sample per patient). Block IDs can be included as
#'   either fixed effects or random intercept terms; this choice is specified using the
#'   argument \code{block_IDs_type}. Note that some testing methods (e.g.
#'   'diffcyt-DA-limma' and 'diffcyt-DS-med') require random intercept block IDs to be
#'   provided directly instead.
#' 
#' @param block_IDs_type (Optional) Whether block IDs should be included in the model
#'   formula as fixed effects or random intercept terms. Options are 'fixed' and 'random'.
#'   Default = 'fixed'.
#' 
#' @param sample_IDs (Optional) Vector or factor of sample IDs. Sample IDs are included as
#'   random intercepts in the model formula, to account for overdispersion typically seen
#'   in high-dimensional cytometry data. This is also known as an 'observation-level
#'   random effect'; for more details see Nowicka et al. (2017), \emph{F1000Research}.
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
                          batch_IDs = NULL, covariates = NULL, 
                          block_IDs = NULL, block_IDs_type = c("fixed", "random"), 
                          sample_IDs = NULL) {
  
  if (!is.factor(group_IDs)) {
    group_IDs <- factor(group_IDs, levels = unique(group_IDs))
  }
  if (!is.null(batch_IDs) & !is.factor(batch_IDs)) {
    batch_IDs <- factor(batch_IDs, levels = unique(batch_IDs))
  }
  if (!is.null(covariates) & !(is.matrix(covariates) & is.numeric(covariates))) {
    stop("'covariates' must be provided as a numeric matrix, with one column for each covariate")
  }
  if (!is.null(block_IDs) & !is.factor(block_IDs)) {
    block_IDs <- factor(block_IDs, levels = unique(block_IDs))
  }
  if (!is.null(sample_IDs) & !is.factor(sample_IDs)) {
    sample_IDs <- factor(sample_IDs, levels = unique(sample_IDs))
  }
  
  # create formula
  
  # fixed effects
  terms <- c("group_IDs", "batch_IDs", "covariates", "block_IDs")
  nulls <- c(is.null(group_IDs), 
             is.null(batch_IDs), is.null(covariates), 
             is.null(block_IDs) | block_IDs_type == "random")
  
  formula_chr <- paste("y ~", paste(terms[!nulls], collapse = " + "))
  
  # random effects
  if (!is.null(block_IDs) & block_IDs_type == "random") {
    formula_chr <- paste0(formula_chr, " + (1 | block_IDs)")
  }
  if (!is.null(sample_IDs)) {
    formula_chr <- paste0(formula_chr, " + (1 | sample_IDs)")
  }
  
  formula <- as.formula(formula_chr)
  
  # create data frame of variables
  data <- data.frame(group_IDs)
  
  if (!is.null(batch_IDs)) data <- cbind(data, batch_IDs)
  if (!is.null(covariates)) data <- cbind(data, covariates)
  if (!is.null(block_IDs)) data <- cbind(data, block_IDs)
  if (!is.null(sample_IDs)) data <- cbind(data, sample_IDs)
  
  # return as list
  list(formula = formula, data = data)
}


