#' Calculate cluster medians and frequencies
#' 
#' Calculate cluster medians (median functional marker expression) and frequencies (number
#' of cells) by cluster and sample
#' 
#' Calculate median expression of functional markers (cluster medians) and number of cells 
#' (cluster frequencies) by cluster and sample (i.e. for each cluster in each sample).
#' 
#' The cluster frequencies are used as weights in the subsequent statistical tests. The
#' cluster medians are either used for directly testing differences in medians
#' ('diffcyt-med'), or for visualization of results from the other methodologies
#' ('diffcyt-FDA' and 'diffcyt-KS').
#' 
#' Results are returned as a new \code{\link[SummarizedExperiment]{SummarizedExperiment}} 
#' object, where rows = clusters, columns = samples, sheets ('assay' slots) = functional 
#' markers. The additional last sheet ('assay' slot) contains the cluster frequencies.
#' 
#' 
#' @param d_se Transformed data object from previous steps, in 
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}} format, with cluster labels 
#'   added in row meta-data using \code{\link{generateClusters}}.
#' 
#' 
#' @return \code{\link[SummarizedExperiment]{SummarizedExperiment}} object, where rows = 
#'   clusters, columns = samples, sheets ('assay' slots) = functional markers. The 
#'   additional last sheet ('assay' slot) contains the cluster frequencies.
#' 
#' 
#' @importFrom SummarizedExperiment SummarizedExperiment assays rowData colData
#' @importFrom dplyr group_by tally summarize
#' @importFrom tidyr complete
#' @importFrom reshape2 acast
#' @importFrom magrittr '%>%'
#' @importFrom stats median
#' @importFrom methods is
#' 
#' @export
#' 
#' @seealso \code{\link{testDA}} \code{\link{testDE_med}} \code{\link{testDE_FDA}}
#'
#' @examples
#' # See full examples in testing functions: testDA, testDE_med, testDE_FDA
#' 
calcMediansAndFreq <- function(d_se) {
  
  if (!is(d_se, "SummarizedExperiment")) {
    stop("Data object must be a 'SummarizedExperiment'")
  }
  
  if (!("cluster" %in% (colnames(rowData(d_se))))) {
    stop(paste0("Data object does not contain cluster labels. Run 'generateClusters()' ", 
                "to generate cluster labels."))
  }
  
  # calculate cluster frequencies
  
  rowdata_df <- as.data.frame(rowData(d_se))
  
  rowdata_df %>% 
    group_by(cluster, sample) %>% 
    tally %>% 
    complete(sample) -> 
    n_cells
  
  n_cells <- acast(n_cells, cluster ~ sample, value.var = "n", fill = 0)
  
  n_cells <- list(n_cells = n_cells)
  
  # calculate cluster medians
  
  assaydata_mx <- assays(d_se)[[1]]
  
  medians_func <- vector("list", sum(colData(d_se)$is_functional))
  func_names <- as.character(colData(d_se)$markers[colData(d_se)$is_functional])
  names(medians_func) <- func_names
  
  clus <- rowData(d_se)$cluster
  smp <- rowData(d_se)$sample
  
  for (i in seq_along(medians_func)) {
    assaydata_i <- assaydata_mx[, func_names[i], drop = FALSE]
    assaydata_i <- as.data.frame(assaydata_i)
    assaydata_i <- cbind(assaydata_i, sample = smp, cluster = clus)
    colnames(assaydata_i)[1] <- "value"
    
    assaydata_i %>% 
      group_by(cluster, sample) %>% 
      summarize(median = median(value)) -> 
      med
    
    med <- acast(med, cluster ~ sample, value.var = "median", fill = NA)
    
    medians_func[[i]] <- med
  }
  
  # create new SummarizedExperiment
  
  list_all <- c(medians_func, n_cells)
  
  row_data <- data.frame(cluster = factor(sort(unique(clus)), levels = sort(unique(clus))))
  col_data <- data.frame(sample = factor(unique(smp), levels = unique(smp)))
  
  d_clus <- SummarizedExperiment(list_all, rowData = row_data, colData = col_data)
  
  d_clus
}


