#' Calculate empirical cumulative distribution functions (ECDFs)
#' 
#' Calculate empirical cumulative distribution functions (ECDFs) required for
#' methods "diffcyt-FDA", "diffcyt-KS", and "diffcyt-LM"
#' 
#' Methods "diffcyt-FDA", "diffcyt-KS", and "diffcyt-LM" use the empirical cumulative
#' distribution functions (ECDFs) of functional marker expression profiles (for each
#' cluster-sample combination) to test for differential expression of functional markers
#' between sample groups.
#' 
#' This function calculates the ECDFs for each cluster-sample combination for each 
#' functional marker, and evaluates each of them at a set of equally-spaced points. The 
#' ECDF values can then be used in the tests for differential expression.
#' 
#' The ECDF values are returned in a new
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} object, where rows = clusters,
#' columns = samples, sheets ('assay' slots) = functional markers. Note that each "value"
#' is a vector of evaluated ECDF values.
#' 
#' 
#' @param d_se Data object from previous steps, in 
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}} format, containing cluster 
#'   labels as a column in the row meta-data (from \code{\link{generateClusters}}).
#' 
#' @param resolution Resolution for evaluating ECDFs. The value of each ECDF is calculated
#'   at this number of equally spaced points between the maximum and minimum observed 
#'   marker expression values for a given cluster-sample combination. Default = 30.
#' 
#' 
#' @return \code{\link[SummarizedExperiment]{SummarizedExperiment}} object, where rows = 
#'   clusters, columns = samples, sheets ('assay' slots) = functional markers. Each entry 
#'   is a vector of values (i.e. the evaluated ECDF values) with length \code{resolution}.
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
#'
#' @examples
#' # See full examples in testing functions: testDA, testDE_med, testDE_FDA
#' 
calcECDFs <- function(d_se, resolution = 30) {
  
  if (!is(d_se, "SummarizedExperiment")) {
    stop("Data object must be a 'SummarizedExperiment'")
  }
  
  if (!("cluster" %in% (colnames(rowData(d_se))))) {
    stop(paste0("Data object does not contain cluster labels. Run 'generateClusters()' ", 
                "to generate cluster labels."))
  }
  
  # calculate ECDFs
  
  rowdata_df <- as.data.frame(rowData(d_se))
  
  assaydata_mx <- assays(d_se)[[1]]
  
  ECDFs_func <- vector("list", sum(colData(d_se)$is_functional))
  func_names <- as.character(colData(d_se)$markers[colData(d_se)$is_functional])
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
  
  # calculate ECDFs for each functional marker, each cluster-sample combination
  # [to do: could possibly replace loop with 'summarize_each'; but is already fast]
  for (i in seq_along(ECDFs_func)) {
    assaydata_i <- assaydata_mx[, func_names[i], drop = FALSE]
    assaydata_i <- as.data.frame(assaydata_i)
    assaydata_i <- cbind(assaydata_i, sample = smp, cluster = clus)
    colnames(assaydata_i)[1] <- "value"
    
    assaydata_i %>% 
      group_by(cluster, sample) %>% 
      summarize(ECDF = list(evaluate_ecdf(value, s_vals(value)))) -> 
      ECDF
    
    ECDF <- acast(ECDF, cluster ~ sample, value.var = "ECDF", fill = NA)
    
    ECDFs_func[[i]] <- ECDF
  }
  
  # create new SummarizedExperiment
  
  row_data <- data.frame(cluster = factor(sort(unique(clus)), levels = sort(unique(clus))))
  col_data <- data.frame(sample = factor(unique(smp), levels = unique(smp)))
  
  d_ecdfs <- SummarizedExperiment(ECDFs_func, rowData = row_data, colData = col_data)
  
  d_ecdfs
}


