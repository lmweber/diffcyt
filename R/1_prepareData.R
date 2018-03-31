#' Prepare data
#' 
#' Prepare data into format for \code{diffcyt} pipeline
#' 
#' Functions in the \code{diffcyt} analysis pipeline assume that input data is provided as
#' a \code{\link[SummarizedExperiment]{SummarizedExperiment}} object, which contains a
#' single matrix of expression values, together with row and column meta-data.
#' 
#' This function accepts a \code{\link[flowCore]{flowSet}} or a list of
#' \code{\link[flowCore]{flowFrame}}s, \code{data.frames}, or matrices as input (i.e. one
#' \code{flowFrame} or list item per sample). The function then concatenates the data
#' tables into a single matrix of values, and adds row and column meta-data.
#' 
#' Row meta-data should be provided as a data frame named \code{sample_info}, containing
#' columns of relevant sample information such as sample IDs and group IDs. This must
#' contain at least a column named \code{sample}.
#' 
#' Column meta-data should be provided as a data frame named \code{marker_info},
#' containing the following columns of marker information. The column names must be as
#' shown.
#' 
#' \itemize{
#' \item \code{marker_name}: protein marker names
#' \item \code{is_marker}: logical vector indicating whether each column is a marker
#' \item \code{is_type_marker}: logical vector indicating whether each column is a 'cell
#' type' marker (for clustering and testing for differential abundance of cell
#' populations)
#' \item \code{is_state_marker}: logical vector indicating whether each column is a 'cell
#' state' marker (for testing for differential states within cell populations)
#' }
#' 
#' Optionally, random subsampling can be used to select an equal number of cells from each
#' sample (\code{subsampling = TRUE}). This can be useful when there are large differences
#' in total numbers of cells per sample, since it ensures that samples with relatively
#' large numbers of cells do not dominate the clustering. However, subsampling should
#' generally not be used when rare cell populations are of interest, due to the
#' significant loss of information if cells from the rare population are discarded.
#' 
#' 
#' @param d_input Input data. Must be a \code{\link[flowCore]{flowSet}} or list of
#'   \code{\link[flowCore]{flowFrame}}s, \code{data.frames}, or matrices as input (one
#'   \code{flowFrame} or list item per sample).
#' 
#' @param sample_info Data frame of sample information, for example sample IDs and group
#'   IDs. Must contain a column named \code{sample}.
#' 
#' @param marker_info Data frame of marker information for each column. This should
#'   contain columns named \code{marker_name}, \code{is_marker}, \code{is_type_marker},
#'   and \code{is_state_marker}. The first column must contain marker names or column
#'   names; the remaining columns are logical vectors indicating whether each column in
#'   the input data is (i) a protein marker, (ii) a 'cell type' marker, and (iii) a 'cell
#'   state' marker.
#' 
#' @param subsampling Whether to use random subsampling to select an equal number of cells
#'   from each sample. Default = FALSE.
#' 
#' @param n_sub Number of cells to select from each sample by random subsampling, if
#'   \code{subsampling = TRUE}. Default = number of cells in smallest sample.
#' 
#' @param seed_sub Random seed for subsampling. Set to an integer value to generate
#'   reproducible results. Default = \code{NULL}.
#' 
#' 
#' @return \code{d_se}: Returns data as a
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}} containing a single matrix
#'   of data (expression values) in the \code{assays} slot, together with row meta-data
#'   (sample information) and column meta-data (marker information). The
#'   \code{sample_info} data frame is also stored in the \code{metadata} slot, and can be
#'   accessed with \code{metadata(d_se)$sample_info}.
#' 
#' 
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom flowCore flowSet exprs
#' @importFrom methods is as
#' 
#' @export
#' 
#' @examples
#' # For a full workflow example demonstrating the use of each function in the 'diffcyt'
#' # pipeline, see the package vignette.
#' 
#' # Create some random data (without differential signal)
#' cofactor <- 5
#' set.seed(123)
#' d_input <- list(
#'   sample1 = sinh(matrix(rnorm(20000, mean = 0, sd = 1), ncol = 20)) * cofactor, 
#'   sample2 = sinh(matrix(rnorm(20000, mean = 0, sd = 1), ncol = 20)) * cofactor, 
#'   sample3 = sinh(matrix(rnorm(20000, mean = 0, sd = 1), ncol = 20)) * cofactor, 
#'   sample4 = sinh(matrix(rnorm(20000, mean = 0, sd = 1), ncol = 20)) * cofactor
#' )
#' 
#' sample_info <- data.frame(
#'   sample = factor(paste0("sample", 1:4)), 
#'   group = factor(c("group1", "group1", "group2", "group2")), 
#'   stringsAsFactors = FALSE
#' )
#' 
#' marker_info <- data.frame(
#'   marker_name = paste0("marker", 1:20), 
#'   is_marker = rep(TRUE, 20), 
#'   is_type_marker = c(rep(TRUE, 10), rep(FALSE, 10)), 
#'   is_state_marker = c(rep(FALSE, 10), rep(TRUE, 10)), 
#'   stringsAsFactors = FALSE
#' )
#' 
#' # Prepare data
#' d_se <- prepareData(d_input, sample_info, marker_info)
#' 
prepareData <- function(d_input, sample_info, marker_info, 
                        subsampling = FALSE, n_sub = NULL, seed_sub = NULL) {
  
  if (!(is(d_input, "list") | is(d_input, "flowSet"))) {
    stop("Input data must be a 'list' or 'flowSet'")
  }
  
  if (is(d_input, "flowSet")) {
    d_input <- as(d_input, "list")
  }
  
  if (all(sapply(d_input, class) == "flowFrame")) {
    d_ex <- lapply(d_input, exprs)
  } else if (all(sapply(d_input, is.data.frame))) {
    d_ex <- lapply(d_input, as.matrix)
  } else if (all(sapply(d_input, is.matrix))) {
    d_ex <- d_input
  } else {
    stop("input data format not recognized (should be a 'flowSet' or a list of 'flowFrames', data.frames, or matrices)")
  }
  
  if (!(nrow(sample_info) == length(d_ex))) {
    stop("number of rows in 'sample_info' data frame must equal the number of samples")
  }
  
  if (!all(sapply(d_input, function(d) is.null(colnames(d)))) | 
      !all(sapply(d_input, function(d) all(colnames(d) == colnames(d_input[[1]]))))) {
    stop("column (marker) names do not match for all samples")
  }
  
  n_cells <- sapply(d_ex, nrow)
  
  if (subsampling) {
    if (is.null(n_sub)) n_sub <- min(n_cells)
    if (!is.null(seed_sub)) set.seed(seed_sub)
    d_ex <- lapply(d_ex, function(d) d[sample(seq_len(nrow(d)), min(n_sub, nrow(d))), , drop = FALSE])
    n_cells <- sapply(d_ex, nrow)
  }
  
  d_combined <- do.call(rbind, d_ex)
  
  # assume marker names are in first column of 'marker_info'
  colnames(d_combined) <- marker_info[, "marker_name"]
  
  # create row meta-data
  stopifnot(is.data.frame(sample_info))
  
  row_data <- as.data.frame(lapply(sample_info, function(col) {
    as.factor(rep(col, n_cells))
  }))
  
  stopifnot(nrow(row_data) == sum(n_cells))
  
  # create column meta-data
  stopifnot(is.data.frame(marker_info), 
            nrow(marker_info) == ncol(d_combined))
  
  col_data <- marker_info
  
  stopifnot("sample" %in% colnames(sample_info))
  
  # create SummarizedExperiment object
  d_se <- SummarizedExperiment(
    d_combined, 
    rowData = row_data, 
    colData = col_data, 
    metadata = list(sample_info = sample_info)
  )
  
  d_se
}


