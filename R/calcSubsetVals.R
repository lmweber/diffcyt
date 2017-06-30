#' Subset values for each cluster-sample combination
#' 
#' Subset values for each cluster-sample combination for each marker, required for
#' differential expression testing
#' 
#' Subsets the expression values for each cluster-sample combination for each marker, and
#' returns as a new \code{\link{SummarizedExperiment}} object. The values can then be
#' provided to the differential expression testing functions.
#' 
#' The expression values are returned in a new
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} object, where rows = clusters,
#' columns = samples, sheets ('assay' slots) = markers. Note that there is a separate
#' table of values ('assay') for each marker, and each 'value' in the tables consists of a
#' vector of expression values.
#' 
#' 
#' @param d_se Data object from previous steps, in
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}} format, containing cluster
#'   labels as a column in the row meta-data (from \code{\link{generateClusters}}).
#' 
#' 
#' @return \code{\link[SummarizedExperiment]{SummarizedExperiment}} object, where rows =
#'   clusters, columns = samples, sheets ('assay' slots) = markers. Each entry is a vector
#'   of expression values with length equal to the number of cells.
#' 
#' 
#' @importFrom SummarizedExperiment SummarizedExperiment assays rowData colData
#' @importFrom dplyr group_by summarize
#' @importFrom tidyr complete
#' @importFrom reshape2 acast
#' @importFrom magrittr '%>%'
#' @importFrom methods is
#' 
#' @export
#' 
#' @seealso to do
#'
#' @examples
#' # See full examples in testing functions.
#' 
calcSubsetVals <- function(d_se) {
  
  if (!is(d_se, "SummarizedExperiment")) {
    stop("Data object 'd_se' must be a 'SummarizedExperiment'")
  }
  
  if (!("cluster" %in% (colnames(rowData(d_se))))) {
    stop("Data object does not contain cluster labels. Run 'generateClusters' to ", 
         "generate cluster labels.")
  }
  
  # subset data values for each marker
  
  assaydata_mx <- assays(d_se)[[1]]
  
  VALS <- vector("list", sum(colData(d_se)$is_marker_col))
  marker_names <- as.character(colData(d_se)$markers[colData(d_se)$is_marker_col])
  names(VALS) <- marker_names
  
  clus <- rowData(d_se)$cluster
  smp <- rowData(d_se)$sample
  
  for (i in seq_along(VALS)) {
    assaydata_i <- assaydata_mx[, marker_names[i], drop = FALSE]
    assaydata_i <- as.data.frame(assaydata_i)
    assaydata_i <- cbind(assaydata_i, sample = smp, cluster = clus)
    colnames(assaydata_i)[1] <- "value"
    
    assaydata_i %>% 
      group_by(cluster, sample) %>% 
      summarize(vals = list(value)) -> 
      vals
    
    # note: can't fill with NAs or zeros
    vals <- acast(vals, cluster ~ sample, value.var = "vals", fill = NULL)
    
    VALS[[i]] <- vals
  }
  
  # check cluster IDs and sample IDs are identical
  for (i in seq_along(VALS)) {
    if (!all(rownames(VALS[[i]]) == rownames(VALS[[1]]))) {
      stop("Cluster IDs do not match")
    }
    if (!all(colnames(VALS[[i]]) == colnames(VALS[[1]]))) {
      stop("Sample IDs do not match")
    }
  }
  
  # create new SummarizedExperiment (rows = clusters, columns = samples)
  
  row_data <- data.frame(
    cluster = factor(rownames(VALS[[1]]), levels = levels(rowData(d_se)$cluster))
  )
  
  col_data <- data.frame(
    sample = factor(colnames(VALS[[1]]), levels = levels(rowData(d_se)$sample))
  )
  
  d_vals <- SummarizedExperiment(VALS, rowData = row_data, colData = col_data)
  
  d_vals
}


