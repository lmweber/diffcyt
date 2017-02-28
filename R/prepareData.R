#' Prepare data
#' 
#' Prepare data into format required by 'diffcyt' functions
#' 
#' Functions in the 'diffcyt' package assume that input data is provided as a
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}}, which contains a single
#' matrix of data, together with row and column meta-data.
#' 
#' This function accepts a list or \code{\link[flowCore]{flowSet}} as input, concatenates 
#' the data tables into a single matrix, and adds row and column meta-data. Row meta-data
#' contains sample labels and group membership labels; column meta-data contains protein
#' marker names.
#' 
#' 
#' @param d_input Input data. Must be a list or \code{\link[flowCore]{flowSet}}.
#' 
#' @param sample_IDs Vector of sample IDs.
#' 
#' @param group_IDs Vector of group IDs.
#' 
#' 
#' @return d_se Returns data as a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   containing a single matrix of data in the \code{assays} slot, together with row
#'   meta-data (sample and group IDs) and column meta-data (protein marker names).
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
#' # prepare data
#' d_se <- prepareData(d_input, sample_IDs, group_IDs)
prepareData <- function(d_input, sample_IDs, group_IDs) {
  
  if (!(is(d_input, "list") | is(d_input, "flowSet"))) {
    stop("Input data must be a 'list' or 'flowSet'")
  }
  
  if (is(d_input, "flowSet")) {
    d_input <- as(d_input, "list")
  }
  
  d_ex <- lapply(d_input, exprs)
  
  n_cells <- sapply(d_ex, nrow)
  
  d_combined <- do.call(rbind, d_ex)
  
  # create SummarizedExperiment with meta-data
  
  row_data <- data.frame(sample = rep(sample_IDs, n_cells), 
                         group = rep(group_IDs, n_cells))
  row_data$sample <- factor(row_data$sample, levels = unique(row_data$sample))
  row_data$group <- factor(row_data$group, levels = unique(row_data$group))
  
  col_data <- data.frame(marker = colnames(d_combined), row.names = colnames(d_combined))
  
  colnames(d_combined) <- NULL
  
  d_out <- SummarizedExperiment(d_combined, rowData = row_data, colData = col_data)
  
  d_out
}


