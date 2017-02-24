#' Plot abundances for top DA cluster
#' 
#' Basic diagnostic plot showing abundances for top differentially abundant (DA) cluster
#' 
#' Saves a plot showing abundances for the top differentially abundant (DA) cluster 
#' detected by \code{\link{testDA}}.
#' 
#' 
#' @param f1 Output from \code{\link{testDA}}.
#' 
#' @param filename Filename and path to save plot. Default path is current working
#'   directory.
#' 
#' 
#' @importFrom limma topTable
#' @importFrom grDevices pdf
#' @importFrom graphics plot lines
#' 
#' @export
#' 
#' @examples
#' # need to create a small example data set for examples
plotTopDACluster <- function(f1, filename = "./plot_top_DA_cluster.pdf") {
  
  # select top DA cluster
  top_DA_cluster <- rownames(limma::topTable(f1))[1]
  
  # [to do: remove hard-coded parameters and make more generalizable]
  pdf(filename, width = 7, height = 7)
  plot(NA, type = "n", 
       xlim = c(0, 16), ylim = c(0, 1.15), xlab = "sample", ylab = "sqrt(prop)", 
       main = "Top 6 differentially abundant clusters (limma)")
  for (i in 1:6) {
    lines(1:8,  topTable(f1, coef = 2, number = 6)[i, 1:8],  type = "l", col = "red")
    lines(9:16, topTable(f1, coef = 2, number = 6)[i, 9:16], type = "l", col = "blue")
  }
  dev.off()
}


