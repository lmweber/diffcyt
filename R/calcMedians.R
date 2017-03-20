#' Calculate cluster medians
#' 
#' Calculate cluster medians (median functional marker expression per cluster per sample)
#' 
#' Calculate median expression of functional markers (cluster medians) per cluster per
#' sample; i.e. one table of values per functional marker; each containing one value for 
#' each cluster in each sample.
#' 
#' The cluster medians are used in the subsequent statistical tests for method 
#' "diffcyt-med", as well as for visualization of results from all methods.
#' 
#' Results are returned as a new \code{\link[SummarizedExperiment]{SummarizedExperiment}} 
#' object, where rows = clusters, columns = samples, sheets ('assay' slots) = functional 
#' markers.
#' 
#' 
#' @param d_se Data object from previous steps, in 
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}} format, containing cluster 
#'   labels as a column in the row meta-data (from \code{\link{generateClusters}}).
#' 
#' 
#' @return \code{\link[SummarizedExperiment]{SummarizedExperiment}} object, where rows = 
#'   clusters, columns = samples, sheets ('assay' slots) = functional markers.
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
calcMedians <- function(d_se) {
  
  if (!is(d_se, "SummarizedExperiment")) {
    stop("Data object must be a 'SummarizedExperiment'")
  }
  
  if (!("cluster" %in% (colnames(rowData(d_se))))) {
    stop(paste0("Data object does not contain cluster labels. Run 'generateClusters()' ", 
                "to generate cluster labels."))
  }
  
  # calculate cluster medians
  
  assaydata_mx <- assays(d_se)[[1]]
  
  medians_func <- vector("list", sum(colData(d_se)$is_functional))
  func_names <- as.character(colData(d_se)$markers[colData(d_se)$is_functional])
  names(medians_func) <- func_names
  
  clus <- rowData(d_se)$cluster
  smp <- rowData(d_se)$sample
  
  # [to do: could possibly replace loop with 'summarize_each'; but is already fast]
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
  
  row_data <- data.frame(cluster = factor(sort(unique(clus)), levels = sort(unique(clus))))
  col_data <- data.frame(sample = factor(unique(smp), levels = unique(smp)))
  
  d_medians <- SummarizedExperiment(medians_func, rowData = row_data, colData = col_data)
  
  d_medians
}


