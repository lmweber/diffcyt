#' Subset values for each cluster-sample combination
#' 
#' Subset values for each cluster-sample combination, required for differential expression
#' testing method 'diffcyt-KS'
#' 
#' Subsets the expression values for each cluster-sample combination, and returns as a new
#' \code{\link{SummarizedExperiment}} object. The values can then be provided to 
#' \code{\link{testDE_KS}} for differential testing with method 'diffcyt-KS' 
#' (Kolmogorov-Smirnov tests).
#' 
#' The expression values are returned in a new 
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} object, where rows = clusters,
#' columns = samples, sheets ('assay' slots) = markers (non-clustering markers only). Note
#' that there is a separate table of values ('assay') for each marker, and each 'value' in
#' the tables consists of a vector of expression values.
#' 
#' 
#' @param d_se Data object from previous steps, in 
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}} format, containing cluster 
#'   labels as a column in the row meta-data (from \code{\link{generateClusters}}).
#' 
#' 
#' @return \code{\link[SummarizedExperiment]{SummarizedExperiment}} object, where rows = 
#'   clusters, columns = samples, sheets ('assay' slots) = markers (non-clustering markers
#'   only). Each entry is a vector of expression values with length equal to the number of
#'   cells.
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
#' @seealso \code{\link{testDA}} \code{\link{testDE_med}} \code{\link{testDE_FDA}} 
#'   \code{\link{testDE_KS}} \code{\link{testDE_LM}}
#'
#' @examples
#' # See full examples in testing functions: testDA, testDE_med, testDE_FDA, testDE_KS,
#' # testDE_LM
#' 
subsetVals <- function(d_se) {
  
  if (!is(d_se, "SummarizedExperiment")) {
    stop("Data object 'd_se' must be a 'SummarizedExperiment'")
  }
  
  if (!("cluster" %in% (colnames(rowData(d_se))))) {
    stop(paste0("Data object does not contain cluster labels. Run 'generateClusters' to ", 
                "generate cluster labels."))
  }
  
  # subset data values for each functional (non-clustering) marker
  
  assaydata_mx <- assays(d_se)[[1]]
  
  vals_func <- vector("list", sum(colData(d_se)$is_DE_col))
  func_names <- as.character(colData(d_se)$markers[colData(d_se)$is_DE_col])
  names(vals_func) <- func_names
  
  clus <- rowData(d_se)$cluster
  smp <- rowData(d_se)$sample
  
  for (i in seq_along(vals_func)) {
    assaydata_i <- assaydata_mx[, func_names[i], drop = FALSE]
    assaydata_i <- as.data.frame(assaydata_i)
    assaydata_i <- cbind(assaydata_i, sample = smp, cluster = clus)
    colnames(assaydata_i)[1] <- "value"
    
    assaydata_i %>% 
      group_by(cluster, sample) %>% 
      summarize(vals = list(value)) -> 
      vals
    
    # note: can't fill with NAs or zeros
    vals <- acast(vals, cluster ~ sample, value.var = "vals", fill = NULL)
    
    vals_func[[i]] <- vals
  }
  
  # check cluster IDs and sample IDs are identical
  for (i in seq_along(vals_func)) {
    if (!all(rownames(vals_func[[i]]) == rownames(vals_func[[1]]))) {
      stop("Cluster IDs do not match")
    }
    if (!all(colnames(vals_func[[i]]) == colnames(vals_func[[1]]))) {
      stop("Sample IDs do not match")
    }
  }
  
  # create new SummarizedExperiment (rows = clusters, columns = samples)
  
  row_data <- data.frame(cluster = factor(rownames(vals_func[[1]]), 
                                          levels = levels(rowData(d_se)$cluster)))
  col_data <- data.frame(sample = factor(colnames(vals_func[[1]]), 
                                         levels = levels(rowData(d_se)$sample)))
  
  d_vals <- SummarizedExperiment(vals_func, rowData = row_data, colData = col_data)
  
  d_vals
}


