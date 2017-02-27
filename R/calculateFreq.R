#' Calculate cluster frequencies
#' 
#' Calculate cluster frequencies (number of cells per cluster)
#' 
#' Calculate number of cells per cluster (cluster frequencies or abundances). This is 
#' required for subsequent differential tests, e.g. differential abundance of clusters or
#' differential expression of functional markers using frequencies as weights.
#' 
#' This function should be run after generating clusters with 
#' \code{\link{generateClusters}}.
#' 
#' 
#' @param d_transf Transformed input data from the previous steps. This should be in the 
#'   form of a \code{\link[flowCore]{flowSet}} object from the 
#'   \code{\link[flowCore]{flowCore}} package.
#' 
#' @param clus Cluster labels for individual cells, from \code{\link{generateClusters}}.
#' 
#' 
#' @return List containing:
#' \itemize{
#' \item tbl_freq: Matrix of cluster frequencies or abundances (rows = clusters, columns 
#' = samples).
#' \item tbl_prop: Matrix of cluster proportions (rows = clusters, columns = samples).
#' }
#' 
#' 
#' @importFrom stats model.matrix
#' @importFrom methods as is
#' 
#' @export
#' 
#' @examples
#' # need to create a small example data set for examples
calculateFreq <- function(d_transf, clus) {
  
  sample_IDs <- rownames(flowCore::phenoData(d_transf))
  
  # number of cells per sample
  n_cells <- sapply(as(d_transf, "list"), nrow)
  
  stopifnot(all(sample_IDs == gsub("\\.[a-z]+$", "", names(n_cells))))  # [to do: generalize to remove regular expression]
  stopifnot(length(sample_IDs) == length(n_cells))
  stopifnot(length(clus) == sum(n_cells))
  
  # calculate table of frequencies
  samp <- rep(sample_IDs, n_cells)
  stopifnot(length(samp) == length(clus))
  # rows = clusters, columns = samples
  tbl_freq <- table(cluster = clus, sample = samp)
  # rearrange columns ('table()' sorts alphabetically; want original order instead)
  tbl_freq <- tbl_freq[, sample_IDs]
  stopifnot(all(colnames(tbl_freq) == sample_IDs))
  
  # table of proportions
  tbl_prop <- t(t(tbl_freq) / colSums(tbl_freq))  # transpose because vector wraps by column
  stopifnot(all(colSums(tbl_prop) == 1))
  
  list(tbl_freq = tbl_freq, tbl_prop = tbl_prop)
}


