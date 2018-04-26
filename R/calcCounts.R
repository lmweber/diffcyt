#' Calculate cluster cell counts
#' 
#' Calculate number of cells per cluster-sample combination
#' 
#' Calculate number of cells per cluster-sample combination (referred to as cluster cell
#' 'counts', 'abundances', or 'frequencies').
#' 
#' The cluster cell counts are required for testing for differential abundance of cell
#' populations, and are also used for weights and filtering when testing for differential
#' states within cell populations.
#' 
#' Results are returned as a new \code{\link{SummarizedExperiment}} object, where rows =
#' clusters, columns = samples, \code{assay} = values (counts). (Note that this structure
#' differs from the input data object.)
#' 
#' 
#' @param d_se Data object from previous steps, in \code{\link{SummarizedExperiment}}
#'   format, containing cluster labels as a column in the row meta-data (from
#'   \code{\link{generateClusters}}).
#' 
#' 
#' @return \code{d_counts}: \code{\link{SummarizedExperiment}} object, where rows =
#'   clusters, columns = samples, \code{assay} = values (counts).
#' 
#' 
#' @importFrom SummarizedExperiment SummarizedExperiment rowData colData
#' @importFrom dplyr group_by tally summarize
#' @importFrom tidyr complete
#' @importFrom reshape2 acast
#' @importFrom magrittr '%>%'
#' @importFrom methods is
#' 
#' @export
#' 
#' @examples
#' # For a complete workflow example demonstrating each step in the 'diffcyt' pipeline, 
#' # see the package vignette.
#' 
#' # Function to create random data (one sample)
#' d_random <- function(n = 20000, mean = 0, sd = 1, ncol = 20, cofactor = 5) {
#'   d <- sinh(matrix(rnorm(n, mean, sd), ncol = ncol)) * cofactor
#'   colnames(d) <- paste0("marker", sprintf("%02d", 1:ncol))
#'   d
#' }
#' 
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
#'   channel_name = paste0("channel", sprintf("%03d", 1:20)), 
#'   marker_name = paste0("marker", sprintf("%02d", 1:20)), 
#'   marker_class = factor(c(rep("type", 10), rep("state", 10)), 
#'                         levels = c("type", "state", "none")), 
#'   stringsAsFactors = FALSE
#' )
#' 
#' # Prepare data
#' d_se <- prepareData(d_input, experiment_info, marker_info)
#' 
#' # Transform data
#' d_se <- transformData(d_se)
#' 
#' # Generate clusters
#' d_se <- generateClusters(d_se)
#' 
#' # Calculate counts
#' d_counts <- calcCounts(d_se)
#' 
calcCounts <- function(d_se) {
  
  if (!is(d_se, "SummarizedExperiment")) {
    stop("Data object must be a 'SummarizedExperiment'")
  }
  
  if (!("cluster_id" %in% (colnames(rowData(d_se))))) {
    stop("Data object does not contain cluster labels. Run 'generateClusters' to generate cluster labels.")
  }
  
  # calculate cluster cell counts
  
  rowdata_df <- as.data.frame(rowData(d_se))
  
  rowdata_df %>% 
    group_by(cluster_id, sample_id) %>% 
    tally %>% 
    complete(sample_id) -> 
    n_cells
  
  n_cells <- acast(n_cells, cluster_id ~ sample_id, value.var = "n", fill = 0)
  
  # fill in any missing clusters
  if (nrow(n_cells) < nlevels(rowData(d_se)$cluster_id)) {
    ix_missing <- which(!(levels(rowData(d_se)$cluster_id) %in% rownames(n_cells)))
    n_cells_tmp <- matrix(0, nrow = length(ix_missing), ncol = ncol(n_cells))
    rownames(n_cells_tmp) <- ix_missing
    n_cells <- rbind(n_cells, n_cells_tmp)
    # re-order rows
    n_cells <- n_cells[order(as.numeric(rownames(n_cells))), , drop = FALSE]
  }
  
  n_cells_total <- rowSums(n_cells)
  
  # create new SummarizedExperiment (with rows = clusters)
  
  row_data <- data.frame(
    cluster_id = factor(rownames(n_cells), levels = levels(rowData(d_se)$cluster_id)), 
    n_cells = n_cells_total, 
    stringsAsFactors = FALSE
  )
  
  col_data <- metadata(d_se)$experiment_info
  
  # rearrange sample order to match 'experiment_info'
  n_cells <- n_cells[, match(col_data$sample_id, colnames(n_cells)), drop = FALSE]
  stopifnot(all(col_data$sample_id == colnames(n_cells)))
  
  d_counts <- SummarizedExperiment(
    assays = list(counts = n_cells), 
    rowData = row_data, 
    colData = col_data
  )
  
  d_counts
}


