#' Create design matrix
#' 
#' Create design matrix for model fitting
#' 
#' Creates a design matrix for model fitting. Analogous to \code{\link{createFormula}},
#' but returns a design matrix instead of model formula.
#' 
#' The design matrix can then be provided to differential testing functions that require a
#' design matrix, together with the data object and contrast matrix.
#' 
#' The \code{group_IDs} input specifies the groups for differential testing. This can be
#' provided as a vector or factor. The first level of the factor or first entry of the
#' vector is used as the 'reference level' for differential testing. To set a different
#' reference level, re-order the levels using \code{\link[stats]{relevel}} or
#' \code{\link[base]{factor}}.
#' 
#' Note that block IDs (e.g. for paired designs) cannot be provided here, since they would
#' then be included in the models as fixed effects. They can instead be provided to
#' \code{\link{createFormula}} to create a model formula instead of a design matrix, or
#' directly to the model fitting function (for some methods only, e.g.
#' 'diffcyt-DA-limma').
#' 
#' 
#' @param group_IDs Vector or factor of group membership labels for each sample (e.g.
#'   diseased vs. healthy, or treated vs. untreated). Vectors are converted to factors
#'   internally. The first level of the factor (e.g. healthy) or first entry of the vector
#'   will be used as the reference level for differential testing. To re-order factor
#'   levels, use \code{\link[stats]{relevel}} or \code{\link[base]{factor}}.
#' 
#' @param batch_IDs (Optional) Vector or factor of batch IDs. Batch IDs are included as
#'   columns of indicator variables in the design matrix. This allows batch effects to be
#'   estimated during model fitting, and taken into account during inference on the group
#'   ID parameters of interest.
#' 
#' @param covariates (Optional) Numeric matrix of continuous covariates. Covariates are
#'   included as columns in the design matrix. This allows their effects to be estimated
#'   during model fitting, and taken into account during inference on the group ID
#'   parameters of interest.
#' 
#' 
#' @return Returns a design matrix (numeric matrix), with one row per sample, and one
#'   model parameter per column.
#' 
#' 
#' @importFrom stats as.formula model.matrix
#' @importFrom methods is
#' 
#' @export
#' 
#' @seealso to do
#' 
#' @examples
#' # to do
#' 
createDesignMatrix <- function(group_IDs, batch_IDs = NULL, covariates = NULL) {
  
  if (!is.factor(group_IDs)) {
    group_IDs <- factor(group_IDs, levels = unique(group_IDs))
  }
  if (!is.null(batch_IDs) & !is.factor(batch_IDs)) {
    batch_IDs <- factor(batch_IDs, levels = unique(batch_IDs))
  }
  if (!is.null(covariates) & !(is.matrix(covariates) & is.numeric(covariates))) {
    stop("'covariates' must be provided as a numeric matrix, with one column for each covariate")
  }
  
  # create design matrix
  terms <- c("group_IDs", "batch_IDs", "covariates")
  nulls <- c(is.null(group_IDs), is.null(batch_IDs), is.null(covariates))
  
  formula <- as.formula(paste("~", paste(terms[!nulls], collapse = " + ")))
  
  design <- model.matrix(formula)
  
  design
}


