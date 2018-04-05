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
#' Results are returned as a new \linkS4class{SummarizedExperiment} object, where rows =
#' clusters, columns = samples, \code{assay} = values (counts). (Note that this structure
#' differs from the input data object.)
#' 
#' 
#' @param d_se Data object from previous steps, in \linkS4class{SummarizedExperiment}
#'   format, containing cluster labels as a column in the row meta-data (from
#'   \code{\link{generateClusters}}).
#' 
#' 
#' @return \code{d_counts}: \linkS4class{SummarizedExperiment} object, where rows =
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
#' sample_info <- data.frame(
#'   sample = factor(paste0("sample", 1:4)), 
#'   group = factor(c("group1", "group1", "group2", "group2")), 
#'   stringsAsFactors = FALSE
#' )
#' 
#' marker_info <- data.frame(
#'   marker_name = paste0("marker", sprintf("%02d", 1:20)), 
#'   is_marker = rep(TRUE, 20), 
#'   marker_type = factor(c(rep("cell_type", 10), rep("cell_state", 10)), 
#'                        levels = c("cell_type", "cell_state", "none")), 
#'   stringsAsFactors = FALSE
#' )
#' 
#' # Prepare data
#' d_se <- prepareData(d_input, sample_info, marker_info)
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
  
  if (!("cluster" %in% (colnames(rowData(d_se))))) {
    stop("Data object does not contain cluster labels. Run 'generateClusters' to generate cluster labels.")
  }
  
  # calculate cluster cell counts
  
  rowdata_df <- as.data.frame(rowData(d_se))
  
  rowdata_df %>% 
    group_by(cluster, sample) %>% 
    tally %>% 
    complete(sample) -> 
    n_cells
  
  n_cells <- acast(n_cells, cluster ~ sample, value.var = "n", fill = 0)
  
  # fill in any missing clusters
  if (nrow(n_cells) < nlevels(rowData(d_se)$cluster)) {
    ix_missing <- which(!(levels(rowData(d_se)$cluster) %in% rownames(n_cells)))
    n_cells_tmp <- matrix(0, nrow = length(ix_missing), ncol = ncol(n_cells))
    rownames(n_cells_tmp) <- ix_missing
    n_cells <- rbind(n_cells, n_cells_tmp)
    # re-order rows
    n_cells <- n_cells[order(as.numeric(rownames(n_cells))), , drop = FALSE]
  }
  
  n_cells_total <- rowSums(n_cells)
  
  # create new SummarizedExperiment (with rows = clusters)
  
  row_data <- data.frame(
    cluster = factor(rownames(n_cells), levels = levels(rowData(d_se)$cluster)), 
    n_cells = n_cells_total, 
    stringsAsFactors = FALSE
  )
  
  col_data <- metadata(d_se)$sample_info
  
  # rearrange sample order to match 'sample_info'
  n_cells <- n_cells[, match(col_data$sample, colnames(n_cells)), drop = FALSE]
  stopifnot(all(col_data$sample == colnames(n_cells)))
  
  d_counts <- SummarizedExperiment(
    n_cells, 
    rowData = row_data, 
    colData = col_data
  )
  
  d_counts
}


