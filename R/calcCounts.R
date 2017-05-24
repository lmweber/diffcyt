#' Calculate cluster cell counts
#' 
#' Calculate number of cells per cluster per sample (i.e. cell counts / abundances /
#' frequencies per cluster-sample combination)
#' 
#' Calculate number of cells per cluster per sample (referred to as cluster cell 'counts',
#' 'abundances', or 'frequencies').
#' 
#' The cluster cell counts are required for the differential abundance tests, and are also
#' used as weights in the differential expression tests and the calculation of false
#' discovery rates (FDRs).
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
#' @seealso \code{\link{testDA}} \code{\link{testDE_med}} \code{\link{testDE_FDA}}
#'   \code{\link{testDE_KS}} \code{\link{testDE_LM}}
#'
#' @examples
#' # See full examples in testing functions: testDA, testDE_med, testDE_FDA, testDE_KS,
#' # testDE_LM
#' 
calcCounts <- function(d_se) {
  
  if (!is(d_se, "SummarizedExperiment")) {
    stop("Data object must be a 'SummarizedExperiment'")
  }
  
  if (!("cluster" %in% (colnames(rowData(d_se))))) {
    stop(paste0("Data object does not contain cluster labels. Run 'generateClusters' to ", 
                "generate cluster labels."))
  }
  
  # calculate cluster cell counts
  
  rowdata_df <- as.data.frame(rowData(d_se))
  
  rowdata_df %>% 
    group_by(cluster, sample) %>% 
    tally %>% 
    complete(sample) -> 
    n_cells
  
  n_cells <- acast(n_cells, cluster ~ sample, value.var = "n", fill = 0)
  
  n_cells_total <- rowSums(n_cells)
  
  # create new SummarizedExperiment (with rows = clusters)
  
  row_data <- data.frame(cluster = factor(rownames(n_cells), 
                                          levels = levels(rowData(d_se)$cluster)), 
                         n_cells = n_cells_total)
  
  col_data <- data.frame(sample = factor(colnames(n_cells), 
                                         levels = levels(rowData(d_se)$sample)))
  
  d_counts <- SummarizedExperiment(n_cells, rowData = row_data, colData = col_data)
  
  d_counts
}


