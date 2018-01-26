#' Prepare data
#' 
#' Prepare data into format required for \code{diffcyt} pipeline
#' 
#' Functions in the \code{diffcyt} analysis pipeline assume that input data is provided as
#' a \code{\link[SummarizedExperiment]{SummarizedExperiment}} object, which contains a
#' single matrix of expression values, together with row and column meta-data.
#' 
#' This function accepts a list or \code{\link[flowCore]{flowSet}} as input (containing
#' one list item or \code{\link[flowCore]{flowFrame}} per sample), concatenates the data
#' tables into a single matrix, and adds row and column meta-data.
#' 
#' Row meta-data contains sample labels and group membership labels. Column meta-data
#' contains protein marker names, and logical entries indicating whether each column is
#' (i) a marker, (ii) a cell type marker (for clustering and testing for differential
#' abundance), and (iii) a functional state marker (for testing for differential
#' functional states).
#' 
#' Optionally, random subsampling can be used to select an equal number of cells from each
#' sample (\code{subsampling = TRUE}). This can be useful when there are large differences
#' in total numbers of cells per sample, since it ensures that samples with relatively
#' large numbers of cells do not dominate the clustering. However, subsampling should
#' generally not be used when rare cell populations are of interest, due to the
#' significant loss of information if cells from the rare population are discarded.
#' 
#' 
#' @param d_input Input data. Must be a list or \code{\link[flowCore]{flowSet}} (one list
#'   item or \code{\link[flowCore]{flowFrame}} per sample).
#' 
#' @param sample_IDs Vector of sample IDs.
#' 
#' @param group_IDs Vector of group IDs.
#' 
#' @param cols_markers Column indices indicating all protein markers.
#' 
#' @param cols_type Column indices indicating cell type markers, to be used for
#'   clustering and testing for differential abundance.
#' 
#' @param cols_state Column indices indicating functional state markers, to be used for
#'   testing for differential functional states.
#' 
#' @param col_names Optional vector of column names; for example, this is useful if the
#'   column names in the input data contain channel names instead of marker names. Default
#'   = column names of input data.
#' 
#' @param subsampling Whether to use random subsampling to select an equal number of cells
#'   from each sample. Default = FALSE.
#' 
#' @param n_sub Number of cells to select from each sample by random subsampling, if
#'   \code{subsampling = TRUE}. Default = number of cells in smallest sample.
#' 
#' 
#' @return d_se Returns data as a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   containing a single matrix of data (expression values) in the \code{assays} slot,
#'   together with row meta-data (sample IDs, group IDs) and column meta-data (protein
#'   marker names, logical vectors for: all markers, cell type markers, and functional
#'   state markers). The \code{group_IDs} vector is also stored in the \code{metadata}
#'   slot, and can be accessed with \code{metadata(d_se)$group_IDs}.
#' 
#' 
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom flowCore flowSet exprs
#' @importFrom methods is as
#' 
#' @export
#' 
#' @examples
#' # A full workflow example demonstrating the use of each function in the 'diffcyt'
#' # pipeline on an experimental data set is available in the package vignette.
#' 
prepareData <- function(d_input, sample_IDs, group_IDs, 
                        cols_markers = NULL, cols_type = NULL, cols_state = NULL, 
                        col_names = NULL, subsampling = FALSE, n_sub = NULL) {
  
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
  
  if (subsampling) {
    if (is.null(n_sub)) n_sub <- min(n_cells)
    d_ex <- lapply(d_ex, function(d) d[sample(seq_len(nrow(d)), n_sub), ])
    n_cells <- sapply(d_ex, nrow)
  }
  
  d_combined <- do.call(rbind, d_ex)
  
  if (!is.null(col_names)) {
    colnames(d_combined) <- col_names
  }
  
  # create row meta-data
  if (!is.factor(sample_IDs)) {
    sample_IDs <- factor(sample_IDs, levels = unique(sample_IDs))
  }
  if (!is.factor(group_IDs)) {
    group_IDs <- factor(group_IDs, levels = unique(group_IDs))
  }
  
  stopifnot(length(sample_IDs) == length(n_cells), length(group_IDs) == length(n_cells))
  
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
  if (!is.null(cols_type)) {
    is_type_col <- empty
    is_type_col[cols_type] <- TRUE
  }
  if (!is.null(cols_state)) {
    is_state_col <- empty
    is_state_col[cols_state] <- TRUE
  }
  
  col_data <- data.frame(
    markers = colnames(d_combined), 
    is_marker_col = is_marker_col, 
    is_type_col = is_type_col, 
    is_state_col = is_state_col, 
    row.names = colnames(d_combined)
  )
  
  colnames(d_combined) <- NULL
  
  # create SummarizedExperiment object
  d_se <- SummarizedExperiment(
    d_combined, 
    rowData = row_data, 
    colData = col_data, 
    metadata = list(group_IDs = group_IDs)
  )
  
  d_se
}


