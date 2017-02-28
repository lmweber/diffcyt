#' Prepare data
#' 
#' Prepare data into format required by 'diffcyt' functions
#' 
#' Functions in the 'diffcyt' package assume that input data is provided as a
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}}, which contains a single
#' matrix of data, together with row and column meta-data.
#' 
#' This function accepts a list or \code{\link[flowCore]{flowSet}} as input, concatenates 
#' the data tables into a single matrix, and adds row and column meta-data.
#' 
#' Row meta-data contains sample labels and group membership labels. Column meta-data 
#' contains protein marker names, and logical entries indicating whether each column is 
#' (i) a marker, (ii) a lineage marker, and (iii) a functional marker.
#' 
#' 
#' @param d_input Input data. Must be a list or \code{\link[flowCore]{flowSet}}.
#' 
#' @param sample_IDs Vector of sample IDs.
#' 
#' @param group_IDs Vector of group IDs.
#' 
#' @param marker_cols Vector of indices of all markers.
#' 
#' @param lineage_cols Vector of indices of lineage markers.
#' 
#' @param functional_cols Vector of indices of functional markers.
#' 
#' 
#' @return d_se Returns data as a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   containing a single matrix of data in the \code{assays} slot, together with row 
#'   meta-data (sample and group IDs) and column meta-data (protein marker names, logical
#'   vectors for: markers, linage markers, functional markers).
#' 
#' 
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom flowCore flowSet exprs
#' @importFrom methods is as
#' 
#' @export
#'
#' @examples
#' library(flowCore)
#' 
#' # filenames
#' files <- list.files(system.file("extdata", package = "diffcyt"), 
#'                     pattern = "\\.fcs$", full.names = TRUE)
#' files_BCRXL <- files[grep("BCRXL", files)]
#' files_ref <- files[grep("ref", files)]
#' 
#' # load data
#' files_load <- c(files_BCRXL, files_ref)
#' d_input <- lapply(files_load, read.FCS, transformation = FALSE, truncate_max_range = FALSE)
#' 
#' # sample IDs and group IDs
#' sample_IDs <- gsub("\\.fcs$", "", basename(files_load))
#' sample_IDs
#' group_IDs <- gsub("^patient[0-9]_", "", sample_IDs)
#' group_IDs
#' 
#' # indices of all marker columns, lineage markers, and functional markers
#' # (see Table 1 in Bruggner et al. 2014)
#' marker_cols <- c(3:4, 7:9, 11:19, 21:22, 24:26, 28:31, 33)
#' lineage_cols <- c(3:4, 9, 11,12,14, 21, 29, 31, 33)
#' functional_cols <- setdiff(marker_cols, lineage_cols)
#' 
#' # prepare data
#' d_se <- prepareData(d_input, sample_IDs, group_IDs, marker_cols, lineage_cols, functional_cols)
prepareData <- function(d_input, sample_IDs, group_IDs, marker_cols, lineage_cols, functional_cols) {
  
  if (!(is(d_input, "list") | is(d_input, "flowSet"))) {
    stop("Input data must be a 'list' or 'flowSet'")
  }
  
  if (is(d_input, "flowSet")) {
    d_input <- as(d_input, "list")
  }
  
  d_ex <- lapply(d_input, exprs)
  
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


