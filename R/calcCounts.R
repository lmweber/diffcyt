#' Calculate cluster counts (frequencies)
#' 
#' Calculate number of cells per cluster per sample (i.e. cluster counts / frequencies / 
#' abundances)
#' 
#' Calculate number of cells per cluster per sample (referred to as cluster 'counts', 
#' 'frequencies', or 'abundances').
#' 
#' The cluster counts are used as weights in the subsequent statistical tests and 
#' calculation of false discovery rates (FDRs).
#' 
#' Results are returned as a new \code{\link[SummarizedExperiment]{SummarizedExperiment}} 
#' object, where rows = clusters, columns = samples, assay = values (counts).
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
#'
#' @examples
#' # See full examples in testing functions: testDA, testDE_med, testDE_FDA
#' 
calcCounts <- function(d_se) {
  
  if (!is(d_se, "SummarizedExperiment")) {
    stop("Data object must be a 'SummarizedExperiment'")
  }
  
  if (!("cluster" %in% (colnames(rowData(d_se))))) {
    stop(paste0("Data object does not contain cluster labels. Run 'generateClusters()' ", 
                "to generate cluster labels."))
  }
  
  # calculate cluster counts
  
  rowdata_df <- as.data.frame(rowData(d_se))
  
  rowdata_df %>% 
    group_by(cluster, sample) %>% 
    tally %>% 
    complete(sample) -> 
    n_cells
  
  n_cells <- acast(n_cells, cluster ~ sample, value.var = "n", fill = 0)
  
  # create new SummarizedExperiment
  
  clus <- rowData(d_se)$cluster
  smp <- rowData(d_se)$sample
  
  row_data <- data.frame(cluster = factor(sort(unique(clus)), levels = sort(unique(clus))))
  col_data <- data.frame(sample = factor(unique(smp), levels = unique(smp)))
  
  d_counts <- SummarizedExperiment(n_cells, rowData = row_data, colData = col_data)
  
  d_counts
}


