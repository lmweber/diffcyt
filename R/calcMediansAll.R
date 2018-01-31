#' Calculate cluster medians across all samples (for plotting)
#' 
#' Calculate cluster medians (median expression of each marker) for each cluster
#' 
#' Calculate median expression of each marker for each cluster across all samples.
#' 
#' The data object is assumed to contain vectors \code{is_marker},
#' \code{is_celltype_marker}, and \code{is_state_marker} in the column meta-data (see
#' \code{\link{prepareData}}). These indicate the sets of all marker columns, cell type
#' marker columns, and state marker columns. Cluster medians are calculated for all
#' markers.
#' 
#' The cluster medians (across all samples) are required for plotting functions.
#' 
#' Variables \code{id_celltype_markers} and \code{id_state_markers} are saved in the
#' \code{metadata} slot of the output object. These can be used to identify the cell type
#' and state markers in later steps of the 'diffcyt' pipeline.
#' 
#' Results are returned as a new \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#' object, where rows = clusters, columns = markers, assay = values (marker expression
#' values). The \code{metadata} slot also contains variables \code{id_celltype_markers}
#' and \code{id_state_markers}, which can be used to identify the sets of cell type and
#' state markers.
#' 
#' 
#' @param d_se Data object from previous steps, in
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}} format, containing cluster
#'   labels as a column in the row meta-data (from \code{\link{generateClusters}}). Column
#'   meta-data is assumed to contain vectors \code{is_marker}, \code{is_celltype_marker},
#'   and \code{is_state_marker}.
#' 
#' 
#' @return \code{d_medians_all}: \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object, where rows = clusters, columns = markers, assay = values (marker expression
#'   values). The \code{metadata} slot contains variables \code{id_celltype_markers} and
#'   \code{id_state_markers}, which can be accessed with
#'   \code{metadata(d_medians)$id_celltype_markers} and
#'   \code{metadata(d_medians)$id_state_markers}.
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
#' @examples
#' # A full workflow example demonstrating the use of each function in the 'diffcyt'
#' # pipeline on an experimental data set is available in the package vignette.
#' 
calcMediansAll <- function(d_se) {
  
  if (!is(d_se, "SummarizedExperiment")) {
    stop("Data object must be a 'SummarizedExperiment'")
  }
  
  if (!("cluster" %in% (colnames(rowData(d_se))))) {
    stop("Data object does not contain cluster labels. Run 'generateClusters' to generate cluster labels.")
  }
  
  # cell type and state markers
  id_celltype_markers <- colData(d_se)$is_celltype_marker[colData(d_se)$is_marker]
  id_state_markers <- colData(d_se)$is_state_marker[colData(d_se)$is_marker]
  
  # calculate cluster medians
  
  marker_vals <- as.data.frame(assay(d_se))[, colData(d_se)$is_marker]
  rowdata_df <- as.data.frame(rowData(d_se))
  
  stopifnot(nrow(marker_vals) == nrow(rowdata_df))
  
  d_all <- cbind(rowdata_df, marker_vals)
  d_all <- melt(d_all, id.vars = seq_len(ncol(rowdata_df)), variable.name = "marker")
  
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
  
  col_data <- colData(d_se)[colData(d_se)$is_marker, ]
  col_data <- col_data[match(colnames(medians), col_data$marker_names), ]
  
  stopifnot(all(col_data$marker_names == colnames(medians)))
  
  metadata <- list(id_celltype_markers = id_celltype_markers, 
                   id_state_markers = id_state_markers)
  
  d_medians_all <- SummarizedExperiment(
    medians, 
    rowData = row_data, 
    colData = col_data, 
    metadata = metadata
  )
  
  d_medians_all
}


