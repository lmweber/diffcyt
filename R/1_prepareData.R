#' Prepare data
#' 
#' Prepare data into format for \code{diffcyt} pipeline
#' 
#' Functions in the \code{diffcyt} analysis pipeline assume that input data is provided as
#' a \code{\link{SummarizedExperiment}} object, which contains a single matrix of
#' expression values, together with row and column meta-data.
#' 
#' This function accepts a \code{\link{flowSet}} or a list of \code{\link{flowFrame}s},
#' \code{data.frames}, or matrices as input (i.e. one \code{flowFrame} or list item per
#' sample). The function then concatenates the data tables into a single matrix of values,
#' and adds row and column meta-data.
#' 
#' Row meta-data should be provided as a data frame named \code{experiment_info},
#' containing columns of relevant experiment information, such as sample IDs and group
#' IDs (for each sample). This must contain at least a column named \code{sample_id}.
#' 
#' Column meta-data should be provided as a data frame named \code{marker_info},
#' containing the following columns of marker information. The column names must be as
#' shown.
#' 
#' \itemize{
#' \item \code{marker_name}: protein marker names (and column names for any other columns)
#' \item \code{marker_class}: factor indicating the protein marker class for each column
#' of data (usually, entries will be either \code{"type"}, \code{"state"}, or
#' \code{"none"})
#' }
#' 
#' The split into 'cell type' and 'cell state' markers is crucial for the analysis. Cell
#' type markers are used to define cell populations by clustering, and to test for
#' differential abundance of cell populations; while cell state markers are used to test
#' for differential states within cell populations.
#' 
#' The optional argument \code{cols_to_include} allows unnecessary columns (e.g. any
#' columns not containing protein markers) to be discarded.
#' 
#' Optionally, random subsampling can be used to select an equal number of cells from each
#' sample (\code{subsampling = TRUE}). This can be useful when there are large differences
#' in total numbers of cells per sample, since it ensures that samples with relatively
#' large numbers of cells do not dominate the clustering. However, subsampling should
#' generally not be used when rare cell populations are of interest, due to the
#' significant loss of information if cells from the rare population are discarded.
#' 
#' 
#' @param d_input Input data. Must be a \code{\link{flowSet}} or list of
#'   \code{\link{flowFrame}s}, \code{\link{DataFrame}s}, \code{data.frames}, or matrices
#'   as input (one \code{flowFrame} or list item per sample).
#' 
#' @param experiment_info \code{data.frame} or \code{\link{DataFrame}} of experiment
#'   information, for example sample IDs and group IDs. Must contain a column named
#'   \code{sample_id}.
#' 
#' @param marker_info \code{data.frame} or \code{\link{DataFrame}} of marker information
#'   for each column of data. This should contain columns named \code{marker_name} and
#'   \code{marker_class}. The columns contain: (i) marker names (and any other column
#'   names); and (ii) a factor indicating the marker class for each column (with entries
#'   \code{"type"}, \code{"state"}, or \code{"none"}).
#' 
#' @param cols_to_include Logical vector indicating which columns to include from the
#'   input data. Default = all columns.
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
#' @return \code{d_se}: Returns data as a \code{\link{SummarizedExperiment}} containing a
#'   single matrix of data (expression values) in the \code{assays} slot, together with
#'   row meta-data (experiment information) and column meta-data (marker information). The
#'   \code{metadata} slot also contains the \code{experiment_info} data frame, and a
#'   vector \code{n_cells} of the number of cells per sample; these can be accessed with
#'   \code{metadata(d_se)$experiment_info} and \code{metadata(d_se)$n_cells}.
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
#' 
#' experiment_info <- data.frame(
#'   sample_id = factor(paste0("sample", 1:4)), 
#'   group_id = factor(c("group1", "group1", "group2", "group2")), 
#'   stringsAsFactors = FALSE
#' )
#' 
#' marker_info <- data.frame(
#'   marker_name = paste0("marker", sprintf("%02d", 1:20)), 
#'   marker_class = factor(c(rep("type", 10), rep("state", 10)), 
#'                         levels = c("type", "state", "none")), 
#'   stringsAsFactors = FALSE
#' )
#' 
#' # Prepare data
#' d_se <- prepareData(d_input, experiment_info, marker_info)
#' 
prepareData <- function(d_input, experiment_info, marker_info, cols_to_include = NULL, 
                        subsampling = FALSE, n_sub = NULL, seed_sub = NULL) {
  
  if (!(is(d_input, "list") | is(d_input, "flowSet"))) {
    stop("Input data must be a 'list' or 'flowSet'")
  }
  
  if (is(d_input, "flowSet")) {
    d_input <- as(d_input, "list")
  }
  
  if (!(all(sapply(d_input, function(d) all(colnames(d) == colnames(d_input[[1]])))) | 
        all(sapply(d_input, function(d) is.null(colnames(d)))))) {
    stop("column (marker) names do not match for all samples")
  }
  
  if (all(sapply(d_input, class) == "flowFrame")) {
    d_ex <- lapply(d_input, exprs)
  } else if (all(sapply(d_input, is.data.frame)) | 
             all(sapply(d_input, class) == "DataFrame")) {
    d_ex <- lapply(d_input, as.matrix)
  } else if (all(sapply(d_input, is.matrix))) {
    d_ex <- d_input
  } else {
    stop("input data format not recognized (should be a 'flowSet' or a list of 'flowFrames', data.frames, or matrices)")
  }
  
  if (!is.null(cols_to_include)) {
    if (!is.logical(cols_to_include)) {
      stop("'cols_to_include' must be a logical vector")
    }
    if (length(cols_to_include) != ncol(d_ex[[1]])) {
      stop("length of 'cols_to_include' does not match number of columns in input data")
    }
    d_ex <- lapply(d_ex, function(e) e[, cols_to_include])
  }
  
  if (!(nrow(marker_info) == ncol(d_ex[[1]]))) {
    stop("number of rows in 'marker_info' data frame must match number of columns in input data (after subsetting with 'cols_to_include')")
  }
  
  if (!(nrow(experiment_info) == length(d_ex))) {
    stop("number of rows in 'experiment_info' data frame must equal the number of samples")
  }
  
  n_cells <- sapply(d_ex, nrow)
  
  if (subsampling) {
    if (is.null(n_sub)) n_sub <- min(n_cells)
    if (!is.null(seed_sub)) set.seed(seed_sub)
    d_ex <- lapply(d_ex, function(d) d[sample(seq_len(nrow(d)), min(n_sub, nrow(d))), , drop = FALSE])
    n_cells <- sapply(d_ex, nrow)
  }
  
  d_combined <- do.call(rbind, d_ex)
  
  # assume marker names are provided in 'marker_info'
  colnames(d_combined) <- marker_info[, "marker_name"]
  
  # create row meta-data
  stopifnot(class(experiment_info) %in% c("data.frame", "DataFrame"), 
            "sample_id" %in% colnames(experiment_info))
  if (class(experiment_info) == "DataFrame") {
    experiment_info <- as.data.frame(experiment_info)
  }
  
  row_data <- as.data.frame(lapply(experiment_info, function(col) {
    as.factor(rep(col, n_cells))
  }))
  
  stopifnot(nrow(row_data) == sum(n_cells))
  
  # create column meta-data
  stopifnot(class(marker_info) %in% c("data.frame", "DataFrame"), 
            nrow(marker_info) == ncol(d_combined))
  if (class(marker_info) == "DataFrame") {
    marker_info <- as.data.frame(marker_info)
  }
  
  col_data <- marker_info
  
  col_data$marker_name <- as.character(col_data$marker_name)
  
  # create SummarizedExperiment object
  d_se <- SummarizedExperiment(
    assays = list(exprs = d_combined), 
    rowData = row_data, 
    colData = col_data, 
    metadata = list(experiment_info = experiment_info, 
                    n_cells = n_cells)
  )
  
  d_se
}


