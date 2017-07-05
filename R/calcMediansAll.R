#' Calculate cluster medians across all samples (for plotting)
#' 
#' Calculate cluster medians (median expression of each marker) for each cluster
#' 
#' Calculate median expression of each marker for each cluster across all samples.
#' 
#' The data object is assumed to contain a vector \code{is_marker_col} in the column
#' meta-data (see \code{\link{prepareData}}), which indicates whether each column
#' represents a marker. Cluster medians are calculated for these columns.
#' 
#' The cluster medians (across all samples) are required for plotting functions.
#' 
#' Results are returned as a new \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#' object, where rows = clusters, columns = markers, assay = values (marker expression
#' values).
#' 
#' 
#' @param d_se Data object from previous steps, in
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}} format, containing cluster
#'   labels as a column in the row meta-data (from \code{\link{generateClusters}}).
#' 
#' 
#' @return \code{\link[SummarizedExperiment]{SummarizedExperiment}} object, where rows =
#'   clusters, columns = markers, assay = values (marker expression values).
#' 
#' 
#' @importFrom SummarizedExperiment SummarizedExperiment assay rowData colData
#' @importFrom dplyr group_by tally summarize
#' @importFrom tidyr complete
#' @importFrom reshape2 melt acast
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
calcMediansAll <- function(d_se) {
  
  if (!is(d_se, "SummarizedExperiment")) {
    stop("Data object must be a 'SummarizedExperiment'")
  }
  
  if (!("cluster" %in% (colnames(rowData(d_se))))) {
    stop("Data object does not contain cluster labels. Run 'generateClusters' to ", 
         "generate cluster labels.")
  }
  
  # calculate cluster medians
  
  marker_vals <- as.data.frame(assay(d_se))
  rowdata_df <- as.data.frame(rowData(d_se))
  
  # remove non-marker values
  marker_vals[, !colData(d_se)$is_marker_col] <- as.numeric(NA)
  
  stopifnot(nrow(marker_vals) == nrow(rowdata_df))
  
  d_all <- cbind(rowdata_df, marker_vals)
  d_all <- melt(d_all, id.vars = 1:ncol(rowdata_df), variable.name = "marker")
  
  d_all %>% 
    group_by(cluster, marker) %>% 
    summarize(median = median(value)) -> 
    medians
  
  medians <- acast(medians, cluster ~ marker, value.var = "median", fill = NA)
  
  # fill in any missing clusters
  if (nrow(medians) < nlevels(rowData(d_se)$cluster)) {
    ix_missing <- which(!(levels(rowData(d_se)$cluster) %in% rownames(medians)))
    medians_tmp <- matrix(NA, nrow = length(ix_missing), ncol = ncol(medians))
    rownames(medians_tmp) <- ix_missing
    medians <- rbind(medians, medians_tmp)
    # re-order rows
    medians <- medians[order(as.numeric(rownames(medians))), ]
  }
  
  # create new SummarizedExperiment (rows = clusters, columns = markers)
  
  row_data <- data.frame(
    cluster = factor(rownames(medians), levels = levels(rowData(d_se)$cluster))
  )
  
  stopifnot(all(colnames(medians) == colData(d_se)$markers), 
            length(colnames(medians)) == length(colData(d_se)$markers))
  
  col_data <- colData(d_se)
  
  d_medians_all <- SummarizedExperiment(medians, rowData = row_data, colData = col_data)
  
  d_medians_all
}


