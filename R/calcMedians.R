#' Calculate cluster medians
#' 
#' Calculate cluster medians (median expression of each marker for each cluster-sample
#' combination)
#' 
#' Calculate median expression of each marker for each cluster-sample combination.
#' 
#' The data object is assumed to contain a vector \code{is_marker_col} in the column
#' meta-data (see \code{\link{prepareData}}), which indicates whether each column
#' represents a marker. Cluster medians are calculated for these columns.
#' 
#' The cluster medians (by sample) are required for differential expression tests.
#' 
#' Results are returned as a new \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#' object, where rows = clusters, columns = samples, sheets ('assay' slots) = markers.
#' Note that there is a separate table of values ('assay') for each marker. The set of
#' functional markers can be identified using the variable 'id_func_markers' saved in the
#' 'metadata' slot.
#' 
#' 
#' @param d_se Data object from previous steps, in
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}} format, containing cluster
#'   labels as a column in the row meta-data (from \code{\link{generateClusters}}).
#' 
#' 
#' @return \code{\link[SummarizedExperiment]{SummarizedExperiment}} object, where rows =
#'   clusters, columns = samples, sheets ('assay' slots) = markers.
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
#' @seealso to do
#'
#' @examples
#' # See full examples in testing functions.
#' 
calcMedians <- function(d_se) {
  
  if (!is(d_se, "SummarizedExperiment")) {
    stop("Data object must be a 'SummarizedExperiment'")
  }
  
  if (!("cluster" %in% (colnames(rowData(d_se))))) {
    stop("Data object does not contain cluster labels. Run 'generateClusters' to ", 
         "generate cluster labels.")
  }
  
  # calculate cluster medians for each marker
  
  assaydata_mx <- assays(d_se)[[1]]
  
  medians <- vector("list", sum(colData(d_se)$is_marker_col))
  marker_names <- as.character(colData(d_se)$markers[colData(d_se)$is_marker_col])
  names(medians) <- marker_names
  # identify functional markers
  id_func_markers <- colData(d_se)$is_DE_col[colData(d_se)$is_marker_col]
  
  clus <- rowData(d_se)$cluster
  smp <- rowData(d_se)$sample
  
  for (i in seq_along(medians)) {
    assaydata_i <- assaydata_mx[, marker_names[i], drop = FALSE]
    assaydata_i <- as.data.frame(assaydata_i)
    assaydata_i <- cbind(assaydata_i, sample = smp, cluster = clus)
    colnames(assaydata_i)[1] <- "value"
    
    assaydata_i %>% 
      group_by(cluster, sample) %>% 
      summarize(median = median(value)) -> 
      med
    
    med <- acast(med, cluster ~ sample, value.var = "median", fill = NA)
    
    medians[[i]] <- med
  }
  
  # check cluster IDs and sample IDs are identical
  for (i in seq_along(medians)) {
    if (!all(rownames(medians[[i]]) == rownames(medians[[1]]))) {
      stop("Cluster IDs do not match")
    }
    if (!all(colnames(medians[[i]]) == colnames(medians[[1]]))) {
      stop("Sample IDs do not match")
    }
  }
  
  # create new SummarizedExperiment (rows = clusters, columns = samples)
  
  row_data <- data.frame(
    cluster = factor(rownames(medians[[1]]), levels = levels(rowData(d_se)$cluster))
  )
  
  stopifnot(all(colnames(medians[[1]]) == levels(rowData(d_se)$sample)))
  
  col_data <- data.frame(
    sample = factor(colnames(medians[[1]]), levels = levels(rowData(d_se)$sample)), 
    group = metadata(d_se)$group_IDs
  )
  
  metadata <- list(id_func_markers = id_func_markers)
  
  d_medians <- SummarizedExperiment(medians, 
                                    rowData = row_data, 
                                    colData = col_data, 
                                    metadata = metadata)
  
  d_medians
}


