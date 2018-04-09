#' Show top clusters or cluster-marker combinations
#' 
#' Show results for top (most highly significant) clusters or cluster-marker combinations
#' 
#' Summary function to display results for top (most highly significant) detected clusters
#' or cluster-marker combinations.
#' 
#' The differential testing functions return results in the form of p-values and adjusted
#' p-values for each cluster (DA tests) or cluster-marker combination (DS tests), which
#' can be used to rank the clusters or cluster-marker combinations by their evidence for
#' differential abundance or differential states. The p-values and adjusted p-values are
#' stored in the \code{rowData} of the output \code{SummarizedExperiment} object generated
#' by the testing functions.
#' 
#' This summary function displays the \code{rowData} from the \code{SummarizedExperiment}
#' output object as a \code{DataFrame}, with clusters or cluster-marker combinations
#' ordered by the adjusted p-values.
#' 
#' 
#' @param res Output object containing results from one of the differential testing
#'   functions, in \code{SummarizedExperiment} format. Differential test results are
#'   stored in \code{rowData}. See \code{\link{testDA_edgeR}}, \code{\link{testDA_voom}},
#'   \code{\link{testDA_GLMM}}, \code{\link{testDS_limma}}, or \code{\link{testDS_LMM}}.
#' 
#' @param order Whether to order results by adjusted p-value. Default = TRUE.
#' 
#' @param all Whether to display all clusters or cluster-marker combinations (instead of
#'   top \code{top_n}). Default = FALSE.
#' 
#' @param top_n Number of clusters or cluster-marker combinations to display (if \code{all
#'   = FALSE}). Default = 20.
#' 
#' 
#' @return Returns a \code{DataFrame} of results for the \code{top_n} clusters or
#'   cluster-marker combinations, ordered by adjusted p-values.
#' 
#' 
#' @importFrom SummarizedExperiment rowData
#' @importFrom utils head
#' 
#' @export
#' 
#' @examples
#' # See the package vignette for a full workflow example showing both types of 
#' # differential analysis (DA and DS), and demonstrating each function in the 'diffcyt' 
#' # pipeline.
#' 
#' # Function to create random data (one sample)
#' d_random <- function(n = 20000, mean = 0, sd = 1, ncol = 20, cofactor = 5) {
#'   sinh(matrix(rnorm(n, mean, sd), ncol = ncol)) * cofactor
#' }
#' # Create random data (without differential signal)
#' set.seed(123)
#' d_input <- list(
#'   sample1 = d_random(), 
#'   sample2 = d_random(), 
#'   sample3 = d_random(), 
#'   sample4 = d_random()
#' )
#' # Add differential signal (for some cells and markers in one group)
#' ix_rows <- 901:1000
#' ix_cols <- c(6:10, 16:20)
#' d_input[[3]][ix_rows, ix_cols] <- d_random(n = 1000, mean = 3, ncol = 10)
#' d_input[[4]][ix_rows, ix_cols] <- d_random(n = 1000, mean = 3, ncol = 10)
#' 
#' experiment_info <- data.frame(
#'   sample_id = factor(paste0("sample", 1:4)), 
#'   group_id = factor(c("group1", "group1", "group2", "group2")), 
#'   stringsAsFactors = FALSE
#' )
#' 
#' marker_info <- data.frame(
#'   marker_name = paste0("marker", sprintf("%02d", 1:20)), 
#'   marker_class = factor(c(rep("cell_type", 10), rep("cell_state", 10)), 
#'                         levels = c("cell_type", "cell_state", "none")), 
#'   stringsAsFactors = FALSE
#' )
#' 
#' # Create design matrix
#' design <- createDesignMatrix(experiment_info, cols_design = 2)
#' # Create contrast matrix
#' contrast <- createContrast(c(0, 1))
#' 
#' # Test for differential abundance (DA) of clusters (using default method 'diffcyt-DA-edgeR')
#' out_DA <- diffcyt(d_input, experiment_info, marker_info, 
#'                   design = design, contrast = contrast, 
#'                   analysis_type = "DA", method_DA = "diffcyt-DA-edgeR")
#' 
#' # Test for differential states (DS) within clusters (using default method 'diffcyt-DS-limma')
#' out_DS <- diffcyt(d_input, experiment_info, marker_info, 
#'                   design = design, contrast = contrast, 
#'                   analysis_type = "DS", method_DS = "diffcyt-DS-limma", 
#'                   plot = FALSE)
#' 
#' # Display results for top DA clusters
#' topClusters(out_DA$res)
#' 
#' # Display results for top DS cluster-marker combinations
#' topClusters(out_DS$res)
#' 
topClusters <- function(res, order = TRUE, all = FALSE, top_n = 20) {
  
  # identify column of adjusted p-values
  ix_p_adj <- which(colnames(rowData(res)) %in% c("FDR", "adj.P.Val", "p_adj"))
  
  res_df <- rowData(res)
  
  if (order) {
    res_df <- res_df[order(res_df[, ix_p_adj]), , drop = FALSE]
  }
  
  if (all) {
    out <- res_df
  } else {
    out <- head(res_df, top_n)
  }
  
  out
}


