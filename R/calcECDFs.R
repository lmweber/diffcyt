#' Calculate empirical cumulative distribution functions (ECDFs)
#' 
#' Calculate empirical cumulative distribution functions (ECDFs) required for differential
#' expression testing
#' 
#' Empirical cumulative distribution functions (ECDFs) of marker expression profiles are
#' used to represent the signal for each cluster-sample combination. The ECDFs can then be
#' used to test for differential expression of (non-clustering) markers between sample
#' groups.
#' 
#' By using the ECDFs, these methods take the full distributions of the marker expression
#' profiles into account. By contrast, most existing approaches rely on only a single
#' feature value or summary statistic for each marker expression profile (such as the
#' median), which discards a significant amount of extra information contained in the rest
#' of the distributions.
#' 
#' This function calculates the ECDFs for each cluster-sample combination for each marker,
#' and evaluates them at a set of equally-spaced points. The evaluated ECDF values can
#' then be used in the differential expression tests (for the non-clustering markers).
#' 
#' The input data object is assumed to contain a vector \code{is_marker_col} in the column
#' meta-data (see \code{\link{prepareData}}), which indicates whether each column
#' represents a marker. ECDFs are calculated for these columns.
#' 
#' The evaluated ECDF values are returned in a new
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} object, where rows = clusters,
#' columns = samples, sheets ('assay' slots) = markers. Note that there is a separate
#' table of values ('assay') for each marker, and each 'value' in the tables consists of a
#' vector of evaluated ECDF values.
#' 
#' 
#' @param d_se Data object from previous steps, in
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}} format, containing cluster
#'   labels as a column in the row meta-data (from \code{\link{generateClusters}}).
#' 
#' @param resolution Resolution for evaluating ECDFs. The value of each ECDF curve is
#'   calculated at this number of equally spaced points between the minimum and maximum
#'   observed (transformed) marker expression values for each cluster-sample combination.
#'   Default = 30.
#' 
#' 
#' @return \code{\link[SummarizedExperiment]{SummarizedExperiment}} object, where rows =
#'   clusters, columns = samples, sheets ('assay' slots) = markers. Each entry is a vector
#'   of values (i.e. the evaluated ECDF values) with length \code{resolution}.
#' 
#' 
#' @importFrom SummarizedExperiment SummarizedExperiment assays rowData colData
#' @importFrom dplyr group_by summarize
#' @importFrom tidyr complete
#' @importFrom reshape2 acast
#' @importFrom magrittr '%>%'
#' @importFrom stats ecdf
#' @importFrom methods is
#' 
#' @export
#' 
#' @seealso to do
#'
#' @examples
#' # See full examples in testing functions.
#' 
calcECDFs <- function(d_se, resolution = 30) {
  
  if (!is(d_se, "SummarizedExperiment")) {
    stop("Data object must be a 'SummarizedExperiment'")
  }
  
  if (!("cluster" %in% (colnames(rowData(d_se))))) {
    stop("Data object does not contain cluster labels. Run 'generateClusters' to ", 
         "generate cluster labels.")
  }
  
  # calculate ECDFs for each marker
  
  assaydata_mx <- assays(d_se)[[1]]
  
  ECDFs <- vector("list", sum(colData(d_se)$is_marker_col))
  marker_names <- as.character(colData(d_se)$markers[colData(d_se)$is_marker_col])
  names(ECDFs) <- marker_names
  
  clus <- rowData(d_se)$cluster
  smp <- rowData(d_se)$sample
  
  # note: 'ecdf' returns a function; evaluate this function at sequence of values 's'
  evaluate_ecdf <- function(vals, s) {
    e <- ecdf(vals)(s)
    names(e) <- s
    e
  }
  s_vals <- function(vals) {
    seq(min(vals), max(vals), length.out = resolution)
  }
  
  for (i in seq_along(ECDFs)) {
    assaydata_i <- assaydata_mx[, marker_names[i], drop = FALSE]
    assaydata_i <- as.data.frame(assaydata_i)
    assaydata_i <- cbind(assaydata_i, sample = smp, cluster = clus)
    colnames(assaydata_i)[1] <- "value"
    
    assaydata_i %>% 
      group_by(cluster, sample) %>% 
      summarize(ECDF = list(evaluate_ecdf(value, s_vals(value)))) -> 
      ECDF
    
    # note: can't fill with NAs or zeros
    ECDF <- acast(ECDF, cluster ~ sample, value.var = "ECDF", fill = NULL)
    
    ECDFs[[i]] <- ECDF
  }
  
  # check cluster IDs and sample IDs are identical
  for (i in seq_along(ECDFs)) {
    if (!all(rownames(ECDFs[[i]]) == rownames(ECDFs[[1]]))) {
      stop("Cluster IDs do not match")
    }
    if (!all(colnames(ECDFs[[i]]) == colnames(ECDFs[[1]]))) {
      stop("Sample IDs do not match")
    }
  }
  
  # create new SummarizedExperiment (rows = clusters, columns = samples)
  
  row_data <- data.frame(
    cluster = factor(rownames(ECDFs[[1]]), levels = levels(rowData(d_se)$cluster))
  )
  
  stopifnot(all(colnames(ECDFs[[1]]) == levels(rowData(d_se)$sample)))
  
  col_data <- data.frame(
    sample = factor(colnames(ECDFs[[1]]), levels = levels(rowData(d_se)$sample)), 
    group = metadata(d_se)$group_IDs
  )
  
  d_ecdfs <- SummarizedExperiment(ECDFs, rowData = row_data, colData = col_data)
  
  d_ecdfs
}


