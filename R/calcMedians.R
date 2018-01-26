#' Calculate cluster medians
#' 
#' Calculate cluster medians (median expression of each marker for each cluster-sample
#' combination)
#' 
#' Calculate median expression of each marker for each cluster-sample combination.
#' 
#' The data object is assumed to contain vectors \code{is_marker}, \code{is_type_marker},
#' and \code{is_state_marker} in the column meta-data (see \code{\link{prepareData}}).
#' These indicate the sets of all marker columns, cell type marker columns, and functional
#' state marker columns. Cluster medians are calculated for all markers.
#' 
#' The cluster medians (by sample) are required for testing for functional states within
#' populations, and for plotting.
#' 
#' Variables \code{id_type_markers} and \code{id_state_markers} are saved in the
#' \code{metadata} slot of the output object. These can be used to identify the cell type
#' and functional state markers in later steps of the 'diffcyt' pipeline.
#' 
#' Results are returned as a new \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#' object, where rows = clusters, columns = samples, sheets ('assay' slots) = markers.
#' Note that there is a separate table of values ('assay') for each marker. The
#' \code{metadata} slot also contains variables \code{id_type_markers} and
#' \code{id_state_markers}, which can be used to identify the sets of cell type and
#' functional state markers.
#' 
#' 
#' @param d_se Data object from previous steps, in
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}} format, containing cluster
#'   labels as a column in the row meta-data (from \code{\link{generateClusters}}). Column
#'   meta-data is assumed to contain vectors \code{is_marker}, \code{is_type_marker}, and
#'   \code{is_state_marker}.
#' 
#' 
#' @return \code{d_medians}: \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object, where rows = clusters, columns = samples, sheets ('assay' slots) = markers.
#'   The \code{metadata} slot contains variables \code{id_type_markers} and
#'   \code{id_state_markers}, which can be accessed with
#'   \code{metadata(d_medians)$id_type_markers} and
#'   \code{metadata(d_medians)$id_state_markers}.
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
#' @examples
#' # A full workflow example demonstrating the use of each function in the 'diffcyt'
#' # pipeline on an experimental data set is available in the package vignette.
#' 
calcMedians <- function(d_se) {
  
  if (!is(d_se, "SummarizedExperiment")) {
    stop("Data object must be a 'SummarizedExperiment'")
  }
  
  if (!("cluster" %in% (colnames(rowData(d_se))))) {
    stop("Data object does not contain cluster labels. Run 'generateClusters' to generate cluster labels.")
  }
  
  # cell type and functional state markers
  id_type_markers <- colData(d_se)$is_type_marker[colData(d_se)$is_marker]
  id_state_markers <- colData(d_se)$is_state_marker[colData(d_se)$is_marker]
  
  # calculate cluster medians for each marker
  
  assaydata_mx <- assays(d_se)[[1]]
  
  medians <- vector("list", sum(colData(d_se)$is_marker))
  marker_names_sub <- as.character(colData(d_se)$marker_names[colData(d_se)$is_marker])
  names(medians) <- marker_names_sub
  
  clus <- rowData(d_se)$cluster
  smp <- rowData(d_se)$sample_IDs
  
  for (i in seq_along(medians)) {
    assaydata_i <- assaydata_mx[, marker_names_sub[i], drop = FALSE]
    assaydata_i <- as.data.frame(assaydata_i)
    assaydata_i <- cbind(assaydata_i, sample_IDs = smp, cluster = clus)
    colnames(assaydata_i)[1] <- "value"
    
    assaydata_i %>% 
      group_by(cluster, sample_IDs) %>% 
      summarize(median = median(value)) -> 
      med
    
    med <- acast(med, cluster ~ sample_IDs, value.var = "median", fill = NA)
    
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
  
  col_data <- metadata(d_se)$sample_info
  col_data <- col_data[match(colnames(medians[[1]]), metadata(d_se)$sample_info$sample_IDs), ]
  
  stopifnot(all(col_data$sample_IDs == colnames(medians[[1]])))
  
  metadata <- list(id_type_markers = id_type_markers, 
                   id_state_markers = id_state_markers)
  
  d_medians <- SummarizedExperiment(
    medians, 
    rowData = row_data, 
    colData = col_data, 
    metadata = metadata
  )
  
  d_medians
}


