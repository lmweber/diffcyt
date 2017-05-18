#' Calculate cluster medians
#' 
#' Calculate cluster medians (median expression of each non-clustering marker for each 
#' cluster-sample combination)
#' 
#' Calculate median expression of each non-clustering marker for each cluster-sample
#' combination.
#' 
#' The data object is assumed to contain a vector \code{is_DE_col} in the column meta-data
#' (see \code{\link{prepareData}}), which indicates whether each column is a 
#' 'non-clustering' marker to be used for differential expression analysis. Cluster 
#' medians are calculated for these markers only.
#' 
#' The cluster medians are used for the differential expression tests, and for
#' visualization of results.
#' 
#' Results are returned as a new \code{\link[SummarizedExperiment]{SummarizedExperiment}} 
#' object, where rows = clusters, columns = samples, sheets ('assay' slots) = markers 
#' (non-clustering markers only). Note that there is a separate table of values ('assay') 
#' for each marker.
#' 
#' 
#' @param d_se Data object from previous steps, in 
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}} format, containing cluster 
#'   labels as a column in the row meta-data (from \code{\link{generateClusters}}).
#' 
#' 
#' @return \code{\link[SummarizedExperiment]{SummarizedExperiment}} object, where rows = 
#'   clusters, columns = samples, sheets ('assay' slots) = markers (non-clustering markers
#'   only).
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
#'   \code{\link{testDE_KS}} \code{\link{testDE_LM}}
#'
#' @examples
#' # See full examples in testing functions: testDA, testDE_med, testDE_FDA, testDE_KS,
#' # testDE_LM
#' 
calcMedians <- function(d_se) {
  
  if (!is(d_se, "SummarizedExperiment")) {
    stop("Data object must be a 'SummarizedExperiment'")
  }
  
  if (!("cluster" %in% (colnames(rowData(d_se))))) {
    stop(paste0("Data object does not contain cluster labels. Run 'generateClusters' to ", 
                "generate cluster labels."))
  }
  
  # calculate cluster medians for each functional (non-clustering) marker
  
  assaydata_mx <- assays(d_se)[[1]]
  
  medians_func <- vector("list", sum(colData(d_se)$is_DE_col))
  func_names <- as.character(colData(d_se)$markers[colData(d_se)$is_DE_col])
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
  
  # check cluster IDs and sample IDs are identical
  for (i in seq_along(medians_func)) {
    if (!all(rownames(medians_func[[i]]) == rownames(medians_func[[1]]))) {
      stop("Cluster IDs do not match")
    }
    if (!all(colnames(medians_func[[i]]) == colnames(medians_func[[1]]))) {
      stop("Sample IDs do not match")
    }
  }
  
  # create new SummarizedExperiment (rows = clusters, columns = samples)
  
  row_data <- data.frame(cluster = factor(rownames(medians_func[[1]]), 
                                          levels = levels(rowData(d_se)$cluster)))
  col_data <- data.frame(sample = factor(colnames(medians_func[[1]]), 
                                         levels = levels(rowData(d_se)$sample)))
  
  d_medians <- SummarizedExperiment(medians_func, rowData = row_data, colData = col_data)
  
  d_medians
}


