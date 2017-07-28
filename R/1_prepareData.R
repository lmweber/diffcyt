#' Prepare data
#' 
#' Prepare data into format required for \code{diffcyt} pipeline
#' 
#' Functions in the \code{diffcyt} analysis pipeline assume that input data is provided as
#' a \code{\link[SummarizedExperiment]{SummarizedExperiment}} object, which contains a
#' single matrix of data values, together with row and column meta-data.
#' 
#' This function accepts a list or \code{\link[flowCore]{flowSet}} as input (containing
#' one list item or \code{\link[flowCore]{flowFrame}} per sample), concatenates the data
#' tables into a single matrix, and adds row and column meta-data.
#' 
#' Row meta-data contains sample labels (e.g. patient IDs) and group membership labels.
#' Column meta-data contains protein marker names, and logical entries indicating whether
#' each column is (i) a marker, (ii) a marker to be used for clustering, and (iii) a
#' marker to be used for differential expression analysis within clusters.
#' 
#' 
#' @param d_input Input data. Must be a list or \code{\link[flowCore]{flowSet}} (one list
#'   item or \code{\link[flowCore]{flowFrame}} per sample).
#' 
#' @param sample_IDs Vector of sample IDs.
#' 
#' @param group_IDs Vector of group IDs. The group IDs also need to be provided separately
#'   to the differential testing functions; they are included here for plotting.
#' 
#' @param cols_markers Column indices indicating all protein markers.
#' 
#' @param cols_clustering Column indices indicating protein markers to be used for
#'   clustering.
#' 
#' @param cols_DE Column indices indicating protein markers to be used for differential
#'   expression analysis.
#' 
#' 
#' @return d_se Returns data as a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   containing a single matrix of data (expression values) in the \code{assays} slot,
#'   together with row meta-data (sample IDs, group IDs) and column meta-data (protein
#'   marker names, logical vectors for: all markers, markers for clustering, markers for
#'   differential expression analysis). The \code{group_IDs} vector is also stored in the
#'   \code{metadata} slot, and can be accessed with \code{metadata(d_se)$group_IDs}.
#' 
#' 
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom flowCore flowSet exprs
#' @importFrom methods is as
#' 
#' @export
#' 
#' @seealso to do
#'
#' @examples
#' # See full examples in testing functions.
#' 
prepareData <- function(d_input, sample_IDs, group_IDs, 
                        cols_markers = NULL, cols_clustering = NULL, cols_DE = NULL) {
  
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
  
  n_cells <- sapply(d_ex, nrow)
  
  d_combined <- do.call(rbind, d_ex)
  
  # create row meta-data
  if (!is.factor(sample_IDs)) {
    sample_IDs <- factor(sample_IDs, levels = unique(sample_IDs))
  }
  if (!is.factor(group_IDs)) {
    group_IDs <- factor(group_IDs, levels = unique(group_IDs))
  }
  
  row_data <- data.frame(
    sample = rep(sample_IDs, n_cells), 
    group = rep(group_IDs, n_cells)
  )
  
  # create column meta-data
  empty <- rep(FALSE, ncol(d_combined))
  if (!is.null(cols_markers)) {
    is_marker_col <- empty
    is_marker_col[cols_markers] <- TRUE
  }
  if (!is.null(cols_clustering)) {
    is_clustering_col <- empty
    is_clustering_col[cols_clustering] <- TRUE
  }
  if (!is.null(cols_DE)) {
    is_DE_col <- empty
    is_DE_col[cols_DE] <- TRUE
  }
  
  col_data <- data.frame(
    markers = colnames(d_combined), 
    is_marker_col = is_marker_col, 
    is_clustering_col = is_clustering_col, 
    is_DE_col = is_DE_col, 
    row.names = colnames(d_combined)
  )
  
  colnames(d_combined) <- NULL
  
  # create SummarizedExperiment object
  d_out <- SummarizedExperiment(d_combined, 
                                rowData = row_data, 
                                colData = col_data, 
                                metadata = list(group_IDs = group_IDs))
  
  d_out
}


