#' Calculate cluster cell counts
#' 
#' Calculate number of cells per cluster per sample
#' 
#' Calculate number of cells per cluster per sample (referred to as cluster cell 'counts',
#' 'abundances', or 'frequencies').
#' 
#' The cluster cell counts are required for testing for differential abundance, and are
#' also used for filtering to improve statistical power when testing for differential
#' functional states.
#' 
#' Results are returned as a new \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#' object, where rows = clusters, columns = samples, assay = values (counts). (Note that
#' this structure differs from the input data object, where rows = cells, and columns =
#' markers.)
#' 
#' 
#' @param d_se Data object from previous steps, in
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}} format, containing cluster
#'   labels as a column in the row meta-data (from \code{\link{generateClusters}}).
#' 
#' 
#' @return \code{\link[SummarizedExperiment]{SummarizedExperiment}} object, where rows =
#'   clusters, columns = samples, assay = values (counts).
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
#' # A full workflow example demonstrating the use of each function in the 'diffcyt'
#' # pipeline on an experimental data set is available in the package vignette.
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
    n_cells <- n_cells[order(as.numeric(rownames(n_cells))), ]
  }
  
  n_cells_total <- rowSums(n_cells)
  
  # create new SummarizedExperiment (with rows = clusters)
  
  row_data <- data.frame(
    cluster = factor(rownames(n_cells), levels = levels(rowData(d_se)$cluster)), 
    n_cells = n_cells_total
  )
  
  stopifnot(all(colnames(n_cells) == levels(rowData(d_se)$sample)))
  
  col_data <- data.frame(
    sample_IDs = levels(rowData(d_se)$sample), 
    group_IDs = metadata(d_se)$group_IDs
  )
  
  d_counts <- SummarizedExperiment(
    n_cells, 
    rowData = row_data, 
    colData = col_data
  )
  
  d_counts
}


