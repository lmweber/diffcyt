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
#' @seealso \code{\link{testDA}} \code{\link{testDE_med}} \code{\link{testDE_FDA}}
#'
#' @examples
#' # See full examples in testing functions: testDA, testDE_med, testDE_FDA
#' 
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


