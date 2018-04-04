#' Calculate cluster medians (across all samples)
#' 
#' Calculate cluster medians across all samples (median expression of each marker for each
#' cluster-marker combination)
#' 
#' Calculate median expression of each marker across all samples, for each cluster-marker
#' combination.
#' 
#' The data object is assumed to contain vectors \code{is_marker} and \code{marker_type}
#' in the column meta-data (see \code{\link{prepareData}}). These indicate (i) whether
#' each column contains a protein marker, and (ii) the protein marker types
#' (\code{"cell_type"}, \code{"cell_state"}, or \code{NA}). Cluster medians are calculated
#' for all markers.
#' 
#' The cluster medians (across all samples) are required for plotting functions.
#' 
#' Variables \code{id_type_markers} and \code{id_state_markers} are saved in the
#' \code{metadata} slot of the output object. These can be used to identify the 'cell
#' type' and 'cell state' markers in the list of \code{assays} in the output
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} object, which is useful in
#' later steps of the 'diffcyt' pipeline.
#' 
#' Results are returned as a new \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#' object, where rows = clusters, columns = markers, assay = values (marker expression
#' values). The \code{metadata} slot also contains variables \code{id_type_markers} and
#' \code{id_state_markers}, which can be used to identify the sets of cell type and cell
#' state markers in the list of \code{assays}.
#' 
#' 
#' @param d_se Data object from previous steps, in
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}} format, containing cluster
#'   labels as a column in the row meta-data (from \code{\link{generateClusters}}). Column
#'   meta-data is assumed to contain vectors \code{is_marker} and \code{marker_type}.
#' 
#' 
#' @return \code{d_medians_all_samples}:
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}} object, where rows =
#'   clusters, columns = markers, assay = values (marker expression values). The
#'   \code{metadata} slot contains variables \code{id_type_markers} and
#'   \code{id_state_markers}, which can be accessed with
#'   \code{metadata(d_medians)$id_type_markers} and
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
#' # For a full workflow example demonstrating the use of each function in the 'diffcyt'
#' # pipeline, see the package vignette.
#' 
#' # Function to create random data (one sample)
#' d_random <- function(n = 20000, mean = 0, sd = 1, ncol = 20, cofactor = 5) {
#'   sinh(matrix(rnorm(n, mean, sd), ncol = ncol)) * cofactor
#' }
#' # Create random data (without differential signal)
#' set.seed(123)
#' d_input <- list(
#'   sample1 = d_random(), 
#'   sample2 = d_random(), 
#'   sample3 = d_random(), 
#'   sample4 = d_random()
#' )
#' 
#' sample_info <- data.frame(
#'   sample = factor(paste0("sample", 1:4)), 
#'   group = factor(c("group1", "group1", "group2", "group2")), 
#'   stringsAsFactors = FALSE
#' )
#' 
#' marker_info <- data.frame(
#'   marker_name = paste0("marker", sprintf("%02d", 1:20)), 
#'   is_marker = rep(TRUE, 20), 
#'   marker_type = factor(c(rep("cell_type", 10), rep("cell_state", 10)), 
#'                        levels = c("cell_type", "cell_state")), 
#'   stringsAsFactors = FALSE
#' )
#' 
#' # Prepare data
#' d_se <- prepareData(d_input, sample_info, marker_info)
#' 
#' # Transform data
#' d_se <- transformData(d_se)
#' 
#' # Generate clusters
#' d_se <- generateClusters(d_se)
#' 
#' # Calculate medians (across all samples)
#' d_medians_all_samples <- calcMediansAllSamples(d_se)
#' 
calcMediansAllSamples <- function(d_se) {
  
  if (!is(d_se, "SummarizedExperiment")) {
    stop("Data object must be a 'SummarizedExperiment'")
  }
  
  if (!("cluster" %in% (colnames(rowData(d_se))))) {
    stop("Data object does not contain cluster labels. Run 'generateClusters' to generate cluster labels.")
  }
  
  # identify 'cell type' and 'cell state' markers in final list of assays
  id_type_markers <- (colData(d_se)$marker_type == "cell_type")[colData(d_se)$is_marker]
  id_state_markers <- (colData(d_se)$marker_type == "cell_state")[colData(d_se)$is_marker]
  
  # calculate cluster medians
  
  marker_vals <- as.data.frame(assay(d_se))[, colData(d_se)$is_marker, drop = FALSE]
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
    medians <- medians[order(as.numeric(rownames(medians))), , drop = FALSE]
  }
  
  # create new SummarizedExperiment (rows = clusters, columns = markers)
  
  row_data <- data.frame(
    cluster = factor(rownames(medians), levels = levels(rowData(d_se)$cluster)), 
    stringAsFactors = FALSE
  )
  
  col_data <- colData(d_se)[colData(d_se)$is_marker, , drop = FALSE]
  
  # rearrange marker order to match 'marker_info'
  medians <- medians[, match(col_data$marker_name, colnames(medians)), drop = FALSE]
  stopifnot(all(col_data$marker_name == colnames(medians)))
  
  metadata <- list(id_type_markers = id_type_markers, 
                   id_state_markers = id_state_markers)
  
  d_medians_all_samples <- SummarizedExperiment(
    medians, 
    rowData = row_data, 
    colData = col_data, 
    metadata = metadata
  )
  
  d_medians_all_samples
}


