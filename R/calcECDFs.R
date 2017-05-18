#' Calculate empirical cumulative distribution functions (ECDFs)
#' 
#' Calculate empirical cumulative distribution functions (ECDFs) required for differential
#' expression testing methods 'diffcyt-FDA', 'diffcyt-KS', and 'diffcyt-LM'
#' 
#' Methods 'diffcyt-FDA', 'diffcyt-KS', and 'diffcyt-LM' use empirical cumulative
#' distribution functions (ECDFs) of the marker expression profiles (for each
#' cluster-sample combination) to test for differential expression of markers between
#' sample groups.
#' 
#' By using the ECDFs, these methods take the full distributions of the marker expression 
#' profiles into account. By contrast, most existing approaches (as well as 'diffcyt-med')
#' rely on only a single feature value or summary statistic for each marker expression
#' profile (e.g. the median), which, discards a significant amount of extra information
#' contained in the rest of the distributions.
#' 
#' This function calculates the ECDFs for each cluster-sample combination for each 
#' non-clustering marker, and evaluates them at a set of equally-spaced points. The 
#' evaluated ECDF values can then be used in the differential expression tests.
#' 
#' The evaluated ECDF values are returned in a new 
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} object, where rows = clusters,
#' columns = samples, sheets ('assay' slots) = markers (non-clustering markers only). Note
#' that there is a separate table of values ('assay') for each marker, and each 'value' in
#' the tables consists of a vector of evaluated ECDF values.
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
#'   clusters, columns = samples, sheets ('assay' slots) = markers (non-clustering markers
#'   only). Each entry is a vector of values (i.e. the evaluated ECDF values) with length 
#'   \code{resolution}.
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
#' @seealso \code{\link{testDA}} \code{\link{testDE_med}} \code{\link{testDE_FDA}}
#'   \code{\link{testDE_KS}} \code{\link{testDE_LM}}
#'
#' @examples
#' # See full examples in testing functions: testDA, testDE_med, testDE_FDA, testDE_KS,
#' # testDE_LM
#' 
calcECDFs <- function(d_se, resolution = 30) {
  
  if (!is(d_se, "SummarizedExperiment")) {
    stop("Data object must be a 'SummarizedExperiment'")
  }
  
  if (!("cluster" %in% (colnames(rowData(d_se))))) {
    stop(paste0("Data object does not contain cluster labels. Run 'generateClusters' to ", 
                "generate cluster labels."))
  }
  
  # calculate ECDFs for each functional (non-clustering) marker
  
  assaydata_mx <- assays(d_se)[[1]]
  
  ECDFs_func <- vector("list", sum(colData(d_se)$is_DE_col))
  func_names <- as.character(colData(d_se)$markers[colData(d_se)$is_DE_col])
  names(ECDFs_func) <- func_names
  
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
  
  for (i in seq_along(ECDFs_func)) {
    assaydata_i <- assaydata_mx[, func_names[i], drop = FALSE]
    assaydata_i <- as.data.frame(assaydata_i)
    assaydata_i <- cbind(assaydata_i, sample = smp, cluster = clus)
    colnames(assaydata_i)[1] <- "value"
    
    assaydata_i %>% 
      group_by(cluster, sample) %>% 
      summarize(ECDF = list(evaluate_ecdf(value, s_vals(value)))) -> 
      ECDF
    
    # note: can't fill with NAs or zeros
    ECDF <- acast(ECDF, cluster ~ sample, value.var = "ECDF", fill = NULL)
    
    ECDFs_func[[i]] <- ECDF
  }
  
  # check cluster IDs and sample IDs are identical
  for (i in seq_along(ECDFs_func)) {
    if (!all(rownames(ECDFs_func[[i]]) == rownames(ECDFs_func[[1]]))) {
      stop("Cluster IDs do not match")
    }
    if (!all(colnames(ECDFs_func[[i]]) == colnames(ECDFs_func[[1]]))) {
      stop("Sample IDs do not match")
    }
  }
  
  # create new SummarizedExperiment (rows = clusters, columns = samples)
  
  row_data <- data.frame(cluster = factor(rownames(ECDFs_func[[1]]), 
                                          levels = levels(rowData(d_se)$cluster)))
  col_data <- data.frame(sample = factor(colnames(ECDFs_func[[1]]), 
                                         levels = levels(rowData(d_se)$sample)))
  
  d_ecdfs <- SummarizedExperiment(ECDFs_func, rowData = row_data, colData = col_data)
  
  d_ecdfs
}


