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
#' @importFrom SummarizedExperiment SummarizedExperiment rowData colData assays
#'   'rowData<-' 'colData<-' 'assays<-'
#' @importFrom stats model.matrix
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
#' # run clustering
#' d_se <- generateClusters(d_se, xdim = 20, ydim = 20, seed = 123, plot = FALSE)
calcMediansAndFreq <- function(d_se) {
  
  if (!is(d_se, "SummarizedExperiment")) {
    stop("Data object must be a 'SummarizedExperiment'")
  }
  
  if (!("cluster" %in% (colnames(rowData(d_se))))) {
    stop(paste0("Data object does not contain cluster labels. Run 'generateClusters()' ", 
                "to generate cluster labels."))
  }
  
  sample_IDs <- rownames(flowCore::phenoData(d_transf))
  
  #clus_all <- do.call("c", clus)
  clus_all <- clus$clus
  
  # number of cells per sample
  n_cells <- sapply(as(d_transf, "list"), nrow)
  
  stopifnot(all(sample_IDs == gsub("\\.[a-z]+$", "", names(n_cells))))  # [to do: generalize to remove regular expression]
  stopifnot(length(sample_IDs) == length(n_cells))
  stopifnot(length(clus_all) == sum(n_cells))
  
  # calculate table of frequencies
  samp <- rep(sample_IDs, n_cells)
  stopifnot(length(samp) == length(clus_all))
  # rows = clusters, columns = samples
  tbl_freq <- table(cluster = clus_all, sample = samp)
  # rearrange columns ('table()' sorts alphabetically; want original order instead)
  tbl_freq <- tbl_freq[, sample_IDs]
  stopifnot(all(colnames(tbl_freq) == sample_IDs))
  
  # table of proportions
  tbl_prop <- t(t(tbl_freq) / colSums(tbl_freq))  # transpose because vector wraps by column
  stopifnot(all(colSums(tbl_prop) == 1))
  
  list(tbl_freq = tbl_freq, tbl_prop = tbl_prop)  # [to do: try to create a SummarizedExperiment containing this]
}


