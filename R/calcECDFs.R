#' Calculate empirical cumulative distribution functions (ECDFs)
#' 
#' Calculate empirical cumulative distribution functions (ECDFs) for FDA-based methods
#' 
#' The functional data analysis (FDA) based methodology in "diffcyt-FDA" uses empirical
#' cumulative distribution functions (ECDFs) of the functional marker expression profiles 
#' by cluster and sample, to test for differential expression between sample groups.
#' 
#' This function calculates the ECDFs for each cluster-sample combination for each 
#' functional marker, and evaluates each of them at a set of equally-spaced points. These 
#' ECDF values can then be used for differential testing with \code{testDE_FDA}.
#' 
#' The ECDF values are returned in a
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} with the same shape as the
#' output from \code{\link{calcMediansAndFreq}}: rows = clusters, columns = samples,
#' sheets ('assay' slots) = functional markers.
#' 
#' 
#' @param d_se Transformed data object from previous steps, in 
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}} format, with cluster labels 
#'   added in row meta-data using \code{\link{generateClusters}}.
#' 
#' @param resolution Resolution for evaluating ECDFs. The value of each ECDF is calculated
#'   at this number of equally spaced points between the maximum and minimum observed
#'   marker expression values for a given cluster-sample combination. Default = 30.
#' 
#' 
#' @return \code{\link[SummarizedExperiment]{SummarizedExperiment}} object, where rows = 
#'   clusters, columns = samples, sheets ('assay' slots) = functional markers. Each entry
#'   is a list of values (i.e. the evaluated ECDF values) with length \code{resolution}.
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
#' @examples
#' library(flowCore)
#' library(limma)
#' 
#' # filenames
#' files <- list.files(system.file("extdata", package = "diffcyt"), 
#'                     pattern = "\\.fcs$", full.names = TRUE)
#' files_BCRXL <- files[grep("BCRXL", files)]
#' files_ref <- files[grep("ref", files)]
#' 
#' # load data
#' files_load <- c(files_BCRXL, files_ref)
#' d_input <- lapply(files_load, read.FCS, transformation = FALSE, truncate_max_range = FALSE)
#' 
#' # sample IDs and group IDs
#' sample_IDs <- gsub("\\.fcs$", "", basename(files_load))
#' sample_IDs
#' group_IDs <- gsub("^patient[0-9]_", "", sample_IDs)
#' group_IDs
#' 
#' # indices of all marker columns, lineage markers, and functional markers
#' # (see Table 1 in Bruggner et al. 2014)
#' marker_cols <- c(3:4, 7:9, 11:19, 21:22, 24:26, 28:31, 33)
#' lineage_cols <- c(3:4, 9, 11,12,14, 21, 29, 31, 33)
#' functional_cols <- setdiff(marker_cols, lineage_cols)
#' 
#' # prepare data
#' d_se <- prepareData(d_input, sample_IDs, group_IDs, marker_cols, lineage_cols, functional_cols)
#' 
#' # transform data
#' d_se <- transformData(d_se, cofactor = 5)
#' 
#' # generate clusters (small 10x10 SOM grid due to small size of example data set)
#' d_se <- generateClusters(d_se, cols_to_use = lineage_cols, xdim = 10, ydim = 10, 
#'                          seed = 123, plot = FALSE)
#' 
#' # calculate cluster medians and frequencies
#' d_clus <- calcMediansAndFreq(d_se)
#' 
#' # calculate ECDFs
#' d_ecdfs <- calcECDFs(d_se)
calcECDFs <- function(d_se, resolution = 30) {
  
  if (!is(d_se, "SummarizedExperiment")) {
    stop("Data object must be a 'SummarizedExperiment'")
  }
  
  if (!("cluster" %in% (colnames(rowData(d_se))))) {
    stop(paste0("Data object does not contain cluster labels. Run 'generateClusters()' ", 
                "to generate cluster labels."))
  }
  
  rowdata_df <- as.data.frame(rowData(d_se))
  
  assaydata_mx <- assays(d_se)[[1]]
  
  ECDFs_func <- vector("list", sum(colData(d_se)$is_functional))
  func_names <- as.character(colData(d_se)$markers[colData(d_se)$is_functional])
  names(ECDFs_func) <- func_names
  
  clus <- rowData(d_se)$cluster
  smp <- rowData(d_se)$sample
  
  # note: 'ecdf' returns a function; evaluate this function at sequence of values 's'
  evaluate_ecdf <- function(vals, s) {
    ecdf(vals)(s)
  }
  s_vals <- function(vals) {
    seq(min(vals), max(vals), length.out = resolution)
  }
  
  # calculate ECDFs for each functional marker; each cluster-sample combination
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


