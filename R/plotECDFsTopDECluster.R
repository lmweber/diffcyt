#' Plot ECDFs for top DE cluster
#' 
#' Basic diagnostic plot showing ECDFs (one per sample) for top differentially expressed 
#' (DE) cluster
#' 
#' Saves a plot showing ECDFs (one per sample) for the top differentially expressed (DE) 
#' cluster detected.
#' 
#' 
#' @param d_ecdfs ECDFs object from \code{\link{calcECDFs}}.
#' 
#' @param path Path to save plot. Default = current working directory.
#' 
#' @param filename Filename for plot.
#' 
#' 
#' @importFrom SummarizedExperiment assays
#' @importFrom limma topTable
#' @importFrom ggplot2 ggplot stat_ecdf aes theme theme_bw unit element_text ggtitle ylab
#'   ggsave
#' 
#' @export
#' 
#' @seealso \code{\link{testDA}} \code{\link{testDE_med}} \code{\link{testDE_FDA}}
#'
#' @examples
#' # See full examples in testing functions: testDA, testDE_med, testDE_FDA
#' 
plotECDFsTopDECluster <- function(d_ecdfs, path = ".", filename = "plot_ECDFs_top_DE_cluster.pdf") {
  
  # [to do: generalize this function; add argument for results object from DE tests, then
  # select clusters automatically]
  
  # select top DE cluster [to do: generalize]
  top_DE_cluster <- "377"
  
  # select marker to plot [to do: generalize]
  top_marker <- "pS6(Yb172)Dd"
  
  d_plot <- assays(d_ecdfs)[[top_marker]][top_DE_cluster, ]
  
  # [to do: add as argument?]
  group <- gsub("^patient[0-9]_", "", names(d_plot))
  group <- factor(group, levels = unique(group))
  
  # fix missing values (no cells in cluster-sample combination)
  resolution <- max(sapply(d_plot, length))
  d_plot <- lapply(d_plot, function(d) {
    if (length(d) == 0) d <- rep(0, resolution)
    d
  })
  
  d_plot <- as.data.frame(d_plot)
  rownames(d_plot) <- NULL
  d_plot <- data.frame(sample = rep(colnames(d_plot), each = nrow(d_plot)), 
                       group = rep(group, each = nrow(d_plot)), 
                       value = unlist(d_plot))
  
  # plot [to do: generalize]
  p <- ggplot(d_plot, aes(value, color = group, group = sample)) + 
    stat_ecdf(geom = "step") + 
    theme(legend.key.size = unit(0.5, "line"), 
          legend.text = element_text(size = 6)) + 
    theme_bw() + 
    ggtitle(paste0("ECDFs for top DE cluster (marker ", top_marker, ")")) + 
    ylab("Cumulative distribution")
  
  ggsave(file.path(path, filename), p, width = 8, height = 7)
}


