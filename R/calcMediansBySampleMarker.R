#' Calculate medians (by sample and marker)
#' 
#' Calculate medians for each sample-marker combination
#' 
#' Calculate overall median marker expression for each sample (i.e. medians for each
#' sample-marker combination).
#' 
#' The data object is assumed to contain a factor \code{marker_class} in the column
#' meta-data (see \code{\link{prepareData}}), which indicates the protein marker class for
#' each column of data (\code{"type"}, \code{"state"}, or \code{"none"}). Cluster medians
#' are calculated for all markers.
#' 
#' The medians by sample and marker are required for plotting purposes.
#' 
#' Variables \code{id_type_markers} and \code{id_state_markers} are saved in the
#' \code{metadata} slot of the output object. These can be used to identify the 'cell
#' type' and 'cell state' markers in the sequence of markers (columns) in the output
#' object, which is useful in later steps of the 'diffcyt' pipeline.
#' 
#' Results are returned as a new \code{\link{SummarizedExperiment}} object, where rows =
#' samples, columns = markers, \code{assay} = values (marker expression values). The
#' \code{metadata} slot also contains variables \code{id_type_markers} and
#' \code{id_state_markers}, which can be used to identify the sets of cell type and cell
#' state markers in the columns.
#' 
#' 
#' @param d_se Data object from previous steps, in \code{\link{SummarizedExperiment}}
#'   format, containing cluster labels as a column in the row meta-data (from
#'   \code{\link{generateClusters}}). Column meta-data is assumed to contain a factor
#'   \code{marker_class}.
#' 
#' 
#' @return \code{d_medians_by_sample_marker}: \code{\link{SummarizedExperiment}} object,
#'   where rows = samples, columns = markers, \code{assay} = values (marker expression
#'   values). The \code{metadata} slot contains variables \code{id_type_markers} and
#'   \code{id_state_markers}, which can be accessed with
#'   \code{metadata(d_medians)$id_type_markers} and
#'   \code{metadata(d_medians)$id_state_markers}.
#' 
#' 
#' @importFrom SummarizedExperiment SummarizedExperiment assays rowData colData
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
#' # For a complete workflow example demonstrating each step in the 'diffcyt' pipeline, 
#' # see the package vignette.
#' 
#' # Function to create random data (one sample)
#' d_random <- function(n = 20000, mean = 0, sd = 1, ncol = 20, cofactor = 5) {
#'   d <- sinh(matrix(rnorm(n, mean, sd), ncol = ncol)) * cofactor
#'   colnames(d) <- paste0("marker", sprintf("%02d", 1:ncol))
#'   d
#' }
#' 
#' # Create random data (without differential signal)
#' set.seed(123)
#' d_input <- list(
#'   sample1 = d_random(), 
#'   sample2 = d_random(), 
#'   sample3 = d_random(), 
#'   sample4 = d_random()
#' )
#' 
#' experiment_info <- data.frame(
#'   sample_id = factor(paste0("sample", 1:4)), 
#'   group_id = factor(c("group1", "group1", "group2", "group2")), 
#'   stringsAsFactors = FALSE
#' )
#' 
#' marker_info <- data.frame(
#'   channel_name = paste0("channel", sprintf("%03d", 1:20)), 
#'   marker_name = paste0("marker", sprintf("%02d", 1:20)), 
#'   marker_class = factor(c(rep("type", 10), rep("state", 10)), 
#'                         levels = c("type", "state", "none")), 
#'   stringsAsFactors = FALSE
#' )
#' 
#' # Prepare data
#' d_se <- prepareData(d_input, experiment_info, marker_info)
#' 
#' # Transform data
#' d_se <- transformData(d_se)
#' 
#' # Generate clusters
#' d_se <- generateClusters(d_se)
#' 
#' # Calculate medians (by sample and marker)
#' d_medians_by_sample_marker <- calcMediansBySampleMarker(d_se)
#' 
calcMediansBySampleMarker <- function(d_se) {
  
  if (!is(d_se, "SummarizedExperiment")) {
    stop("Data object must be a 'SummarizedExperiment'")
  }
  
  is_marker <- colData(d_se)$marker_class != "none"
  
  # identify 'cell type' and 'cell state' markers in final columns
  id_type_markers <- (colData(d_se)$marker_class == "type")[is_marker]
  id_state_markers <- (colData(d_se)$marker_class == "state")[is_marker]
  
  # calculate medians
  
  marker_vals <- as.matrix(assays(d_se)[["exprs"]])[, is_marker, drop = FALSE]
  rowdata_df <- as.data.frame(rowData(d_se))
  
  stopifnot(nrow(marker_vals) == nrow(rowdata_df))
  
  d_all <- cbind(rowdata_df, marker_vals)
  d_all <- melt(d_all, id.vars = seq_len(ncol(rowdata_df)), variable.name = "marker_id")
  
  d_all %>% 
    group_by(sample_id, marker_id, .drop = FALSE) %>% 
    summarize(median = median(value)) -> 
    medians
  
  medians <- acast(medians, sample_id ~ marker_id, value.var = "median", fill = NA)
  
  # fill in any missing clusters
  if (nrow(medians) < nlevels(rowData(d_se)$sample_id)) {
    ix_missing <- which(!(levels(rowData(d_se)$sample_id) %in% rownames(medians)))
    medians_tmp <- matrix(NA, nrow = length(ix_missing), ncol = ncol(medians))
    rownames(medians_tmp) <- ix_missing
    medians <- rbind(medians, medians_tmp)
    # re-order rows
    medians <- medians[order(as.numeric(rownames(medians))), , drop = FALSE]
  }
  
  # create new SummarizedExperiment (rows = samples, columns = markers)
  
  row_data <- data.frame(
    sample_id = factor(rownames(medians), levels = levels(rowData(d_se)$sample_id)), 
    stringsAsFactors = FALSE
  )
  
  col_data <- colData(d_se)[is_marker, , drop = FALSE]
  
  # rearrange marker order to match 'marker_info'
  medians <- medians[, match(col_data$marker_name, colnames(medians)), drop = FALSE]
  stopifnot(all(col_data$marker_name == colnames(medians)))
  
  metadata <- list(id_type_markers = id_type_markers, 
                   id_state_markers = id_state_markers)
  
  d_medians_by_sample_marker <- SummarizedExperiment(
    assays = list(medians_by_sample_marker = medians), 
    rowData = row_data, 
    colData = col_data, 
    metadata = metadata
  )
  
  d_medians_by_sample_marker
}


