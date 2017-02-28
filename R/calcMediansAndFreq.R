#' Calculate cluster medians and frequencies
#' 
#' Calculate cluster medians (median functional marker expression) and frequencies (number
#' of cells) by cluster and sample
#' 
#' Calculate median expression of functional markers (cluster medians) and number of cells 
#' (cluster frequencies) by cluster and sample (i.e. for each cluster in each sample).
#' 
#' The cluster frequencies are used as weights in the subsequent statistical tests. The
#' cluster medians are either used for directly testing differences in medians
#' ('diffcyt-med'), or for visualization of results from the other methodologies
#' ('diffcyt-FDA' and 'diffcyt-KS').
#' 
#' Results are returned as a new \code{\link[SummarizedExperiment]{SummarizedExperiment}} 
#' object, where rows = clusters, columns = samples, sheets ('assay' slots) = functional 
#' markers. The additional last sheet ('assay' slot) contains the cluster frequencies.
#' 
#' 
#' @param d_se Transformed data object from previous steps, in 
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}} format, with cluster labels 
#'   added in row meta-data using \code{\link{generateClusters}}.
#' 
#' 
#' @return \code{\link[SummarizedExperiment]{SummarizedExperiment}} object, where rows = 
#'   clusters, columns = samples, sheets ('assay' slots) = functional markers. The 
#'   additional last sheet ('assay' slot) contains the cluster frequencies.
#' 
#' 
#' @importFrom SummarizedExperiment SummarizedExperiment rowData colData
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
#' library(flowCore)
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
#' d_se <- generateClusters(d_se, xdim = 10, ydim = 10, seed = 123, plot = FALSE)
#' 
#' # calculate cluster medians and frequencies
#' d_clus <- calcMediansAndFreq(d_se)
calcMediansAndFreq <- function(d_se) {
  
  if (!is(d_se, "SummarizedExperiment")) {
    stop("Data object must be a 'SummarizedExperiment'")
  }
  
  if (!("cluster" %in% (colnames(rowData(d_se))))) {
    stop(paste0("Data object does not contain cluster labels. Run 'generateClusters()' ", 
                "to generate cluster labels."))
  }
  
  # calculate cluster frequencies
  
  rowdata_df <- as.data.frame(rowData(d_se))
  
  rowdata_df %>% 
    group_by(cluster, sample) %>% 
    tally %>% 
    complete(sample) -> 
    n_cells
  
  n_cells <- acast(n_cells, cluster ~ sample, value.var = "n", fill = 0)
  
  n_cells <- list(n_cells = n_cells)
  
  # calculate cluster medians
  
  assaydata_mx <- assays(d_se)[[1]]
  
  medians_func <- vector("list", sum(colData(d_se)$is_functional))
  func_names <- as.character(colData(d_se)$markers[colData(d_se)$is_functional])
  names(medians_func) <- func_names
  
  clus <- rowData(d_se)$cluster
  smp <- rowData(d_se)$sample
  
  for (i in seq_along(medians_func)) {
    assaydata_i <- assaydata_mx[, func_names[i], drop = FALSE]
    assaydata_i <- as.data.frame(assaydata_i)
    assaydata_i <- cbind(assaydata_i, sample = smp, cluster = clus)
    colnames(assaydata_i)[1] <- "value"
    
    assaydata_i %>% 
      group_by(cluster, sample) %>% 
      summarize(median = median(value)) -> 
      med
    
    med <- acast(med, cluster ~ sample, value.var = "median", fill = NA)
    
    medians_func[[i]] <- med
  }
  
  # create new SummarizedExperiment
  
  list_all <- c(medians_func, n_cells)
  
  row_data <- data.frame(cluster = factor(sort(unique(clus)), levels = sort(unique(clus))))
  col_data <- data.frame(sample = factor(unique(smp), levels = unique(smp)))
  
  d_clus <- SummarizedExperiment(list_all, rowData = row_data, colData = col_data)
  
  d_clus
}


