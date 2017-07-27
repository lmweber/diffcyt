#' Create model formula and data frame of variables
#' 
#' Create model formula and corresponding data frame of variables for model fitting
#' 
#' Creates a model formula and corresponding data frame of variables for model fitting.
#' Analogous to \code{\link{createDesignMatrix}}, but returns a model formula instead of a
#' design matrix.
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
#'   fixed effects in the model formula. This allows batch effects to be estimated during
#'   model fitting and taken into account during inference on the parameters of interest.
#' 
#' @param covariates (Optional) Numeric matrix of continuous covariates. Covariates are
#'   included as fixed effects in the model formula. This allows their effects to be
#'   estimated during model fitting and taken into account during inference on the
#'   parameters of interest.
#' 
#' @param block_IDs (Optional) Vector or factor of block IDs, e.g. for paired designs
#'   (e.g. one diseased and one healthy sample per patient). Block IDs are included as
#'   random intercepts in the model formula, to improve power during inference on the
#'   parameters of interest.
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
#' @seealso to do
#' 
#' @examples
#' # to do
#' 
createFormula <- function(group_IDs, batch_IDs = NULL, covariates = NULL, 
                          block_IDs = NULL, sample_IDs = NULL) {
  
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
  terms <- c("group_IDs", "batch_IDs", "covariates", "block_IDs", "sample_IDs")
  nulls <- c(is.null(group_IDs), is.null(batch_IDs), is.null(covariates), is.null(block_IDs), is.null(sample_IDs))
  
  formula <- as.formula(paste("y ~", paste(terms[!nulls], collapse = " + ")))
  
  # create data frame of variables
  data <- data.frame(group_IDs)
  
  if (!is.null(batch_IDs)) data <- cbind(data, batch_IDs)
  if (!is.null(covariates)) data <- cbind(data, covariates)
  if (!is.null(block_IDs)) data <- cbind(data, block_IDs)
  if (!is.null(sample_IDs)) data <- cbind(data, sample_IDs)
  
  # return as list
  list(formula = formula, data = data)
}


