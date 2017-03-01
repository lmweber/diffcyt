#' Plot abundances for top DA cluster
#' 
#' Basic diagnostic plot showing abundances for top differentially abundant (DA) cluster
#' 
#' Saves a plot showing abundances for the top differentially abundant (DA) cluster 
#' detected by \code{\link{testDA}}.
#' 
#' 
#' @param res_DA Fitted model objects from differential abundance tests with
#'   \code{\link{testDA}}.
#' 
#' @param path Path to save plot. Default = current working directory.
#' 
#' @param filename Filename for plot.
#' 
#' 
#' @importFrom limma topTable
#' @importFrom grDevices pdf
#' @importFrom graphics plot lines
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
#' # (1) test for differentially abundant (DA) clusters
#' group <- factor(group_IDs, levels = c("ref", "BCRXL"))  # re-level factor to use "ref" as base level
#' res_DA <- testDA(d_clus, group)
#' topTable(res_DA, number = 6)
#' 
#' # plot top differentially abundant (DA) clusters
#' # note there is no evidence for DA in this example data set (data set is too small)
#' # plotTopDAClusters(res_DA)
plotTopDAClusters <- function(res_DA, path = ".", filename = "plot_top_DA_clusters.pdf") {
  
  # [to do: make more flexible]
  pdf(file.path(path, filename), width = 7, height = 7)
  plot(NA, type = "n", 
       xlim = c(1, 16), ylim = c(0, 1.5), xlab = "sample", ylab = "square root proportion", 
       main = "Top 6 differentially abundant clusters (limma)")
  for (i in 1:16) {
    lines(1:8,  topTable(res_DA, coef = 2, number = 6)[i, 1:8],  type = "l", col = "red")
    lines(9:16, topTable(res_DA, coef = 2, number = 6)[i, 9:16], type = "l", col = "blue")
  }
  dev.off()
}


