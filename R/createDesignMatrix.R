#' Create design matrix
#' 
#' Create design matrix for model fitting
#' 
#' Creates a design matrix specifying the models to be fitted. (Alternatively,
#' \code{\link{createFormula}} can be used to generate a model formula instead of a design
#' matrix.)
#' 
#' The design matrix can then be provided to the differential testing functions, together
#' with the data object and contrast matrix.
#' 
#' The \code{group_IDs} input specifies the groups for differential testing. This can be
#' provided as a vector or factor. The first level of the factor or first entry of the
#' vector is used as the 'reference level' for differential testing. To set a different
#' reference level, re-order the levels using \code{\link[stats]{relevel}} or
#' \code{\link[base]{factor}}.
#' 
#' Depending on the experimental design, there may also be block IDs (e.g. patient IDs in
#' a paired design), batch effects, and continuous covariates. These can be included in
#' the design matrix by providing the optional arguments \code{block_IDs},
#' \code{batch_IDs}, and \code{covariates}.
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
#'   designs, such as one diseased and one healthy sample per patient). If provided, block
#'   IDs are included as columns of indicator variables representing fixed effect terms in
#'   the design matrix. This allows their effects to be estimated during model fitting,
#'   and taken into account during inference on the group ID parameters of interest. Note
#'   that block IDs can also be included as random effects (instead of fixed effects) by
#'   using \code{\link{createFormula}}; or, for some methods, by providing them directly
#'   to the differential testing function (\code{\link{testDA_limma}} and
#'   \code{\link{testDS_med}}).
#' 
#' @param batch_IDs (Optional) Vector or factor of batch effect IDs. Batch effect IDs are
#'   included as columns of indicator variables representing fixed effect terms in the
#'   design matrix. This allows batch effects to be estimated during model fitting, and
#'   taken into account during inference on the group ID parameters of interest.
#' 
#' @param covariates (Optional) Numeric matrix of continuous covariates (one column per
#'   covariate). Covariates are included as columns of values in the design matrix. This
#'   allows their effects to be estimated during model fitting, and taken into account
#'   during inference on the group ID parameters of interest.
#' 
#' 
#' @return Returns a design matrix (numeric matrix), with one row per sample, and one
#'   column per model parameter.
#' 
#' 
#' @importFrom stats as.formula model.matrix
#' @importFrom methods is
#' 
#' @export
#' 
#' @examples
#' # A full workflow example demonstrating the use of each function in the 'diffcyt'
#' # pipeline on an experimental data set is available in the package vignette.
#' 
createDesignMatrix <- function(group_IDs, 
                               block_IDs = NULL, batch_IDs = NULL, covariates = NULL) {
  
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
  
  # create design matrix
  terms <- c("group_IDs", "block_IDs", "batch_IDs", "covariates")
  nulls <- c(is.null(group_IDs), is.null(block_IDs), is.null(batch_IDs), is.null(covariates))
  
  formula <- as.formula(paste("~", paste(terms[!nulls], collapse = " + ")))
  
  design <- model.matrix(formula)
  
  design
}


