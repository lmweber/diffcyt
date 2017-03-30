#' Prepare data
#' 
#' Prepare data into format required for \code{diffcyt} pipeline
#' 
#' Functions in the \code{diffcyt} analysis pipeline assume that input data is provided as
#' a \code{\link[SummarizedExperiment]{SummarizedExperiment}}, which contains a single 
#' matrix of data values, together with row and column meta-data.
#' 
#' This function accepts a list or \code{\link[flowCore]{flowSet}} as input (containing
#' one list item or \code{\link[flowCore]{flowFrame}} per sample), concatenates the data
#' tables into a single matrix, and adds row and column meta-data.
#' 
#' Row meta-data contains sample labels and group membership labels. Column meta-data 
#' contains protein marker names, and logical entries indicating whether each column is 
#' (i) a marker, (ii) a lineage marker, and (iii) a functional marker.
#' 
#' 
#' @param d_input Input data. Must be a list or \code{\link[flowCore]{flowSet}} (one list
#'   item or \code{\link[flowCore]{flowFrame}} per sample).
#' 
#' @param sample_IDs Vector of sample IDs.
#' 
#' @param group_IDs Vector of group IDs.
#' 
#' @param marker_cols Vector of indices of all markers.
#' 
#' @param lineage_cols Vector of indices of lineage markers (for clustering).
#' 
#' @param functional_cols Vector of indices of functional markers (for differential
#'   expression analysis).
#' 
#' 
#' @return d_se Returns data as a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   containing a single matrix of data (expression values) in the \code{assays} slot,
#'   together with row meta-data (sample and group IDs) and column meta-data (protein
#'   marker names, logical vectors for: markers, lineage markers, functional markers).
#' 
#' 
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom flowCore flowSet exprs
#' @importFrom methods is as
#' 
#' @export
#' 
#' @seealso \code{\link{testDA}} \code{\link{testDE_med}} \code{\link{testDE_FDA}}
#'
#' @examples
#' # See full examples in testing functions: testDA, testDE_med, testDE_FDA
#' 
prepareData <- function(d_input, sample_IDs, group_IDs, 
                        marker_cols = NULL, lineage_cols = NULL, functional_cols = NULL) {
  
  if (!(is(d_input, "list") | is(d_input, "flowSet"))) {
    stop("Input data must be a 'list' or 'flowSet'")
  }
  
  if (is(d_input, "flowSet")) {
    d_input <- as(d_input, "list")
  }
  
  d_ex <- lapply(d_input, exprs)
  
  if (!(length(sample_IDs) == length(d_ex))) {
    stop("'sample_IDs' vector must have length equal to number of samples")
  }
  if (!(length(group_IDs) == length(d_ex))) {
    stop("'group_IDs' vector must have length equal to number of samples")
  }
  
  n_cells <- sapply(d_ex, nrow)
  
  d_combined <- do.call(rbind, d_ex)
  
  # row meta-data
  row_data <- data.frame(sample = rep(sample_IDs, n_cells), 
                         group = rep(group_IDs, n_cells))
  row_data$sample <- factor(row_data$sample, levels = unique(row_data$sample))
  row_data$group <- factor(row_data$group, levels = unique(row_data$group))
  
  # column meta-data
  is_marker <- is_lineage <- is_functional <- rep(FALSE, ncol(d_combined))
  is_marker[marker_cols] <- TRUE
  is_lineage[lineage_cols] <- TRUE
  is_functional[functional_cols] <- TRUE
  
  col_data <- data.frame(markers = colnames(d_combined), 
                         is_marker = is_marker, 
                         is_lineage = is_lineage, 
                         is_functional = is_functional, 
                         row.names = colnames(d_combined))
  
  colnames(d_combined) <- NULL
  
  # create SummarizedExperiment object
  d_out <- SummarizedExperiment(d_combined, rowData = row_data, colData = col_data)
  
  d_out
}


