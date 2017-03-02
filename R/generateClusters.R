#' Generate clusters
#' 
#' Generate clusters for 'diffcyt' analysis
#' 
#' Runs a clustering algorithm to divide cells into a large number of small clusters, 
#' which will then be further analyzed for differential abundance or differential 
#' expression of functional markers.
#' 
#' Data is assumed to be in the form of a
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}}, and already transformed with
#' \code{\link{transformData}}.
#' 
#' The column meta-data of the input \code{SummarizedExperiment} object is assumed to
#' contain a vector of logical entries (\code{is_lineage}) indicating which columns
#' represent lineage markers. Only the lineage markers will be used for clustering.
#' 
#' We use the \code{\link[FlowSOM]{FlowSOM}} clustering algorithm by default (Van Gassen 
#' et al. 2015, \emph{Cytometry Part A}). We previously showed that \code{FlowSOM} gives
#' very good clustering performance for high-dimensional cytometry data, for both major
#' and rare cell populations, and is extremely fast (Weber and Robinson, 2016, 
#' \emph{Cytometry Part A}).
#' 
#' The clustering is run at high resolution to give a large number of small clusters, to
#' ensure that subtle differential expression signals in small subsets are not missed.
#' 
#' 
#' @param d_se Transformed input data, from \code{\link{transformData}}.
#' 
#' @param cols_to_use Columns to use for clustering. Default = NULL, in which case the
#'   lineage markers from \code{is_lineage} in the column meta-data of \code{d_se} will be
#'   used.
#' 
#' @param xdim Width of grid for self-organizing map for FlowSOM clustering (number of 
#'   clusters = \code{xdim} * \code{ydim}). Default value = 20 (400 clusters).
#' 
#' @param ydim Height of grid for self-organizing map for FlowSOM clustering (number of 
#'   clusters = \code{xdim} * \code{ydim}). Default value = 20 (400 clusters).
#' 
#' @param seed Random seed (set to an integer value for reproducible results). Default 
#'   value = NULL.
#' 
#' @param plot Whether to save plot. Default = TRUE.
#' 
#' @param path Path to save plot.
#' 
#' @param filename Filename for plot.
#' 
#' @param ... Other parameters to pass to FlowSOM (see \code{\link[FlowSOM]{BuildSOM}} 
#'   and \code{\link[FlowSOM]{SOM}} for details.)
#' 
#' 
#' @return Returns the \code{\link[SummarizedExperiment]{SummarizedExperiment}} input 
#'   object with an additional column of row meta-data, containing cluster labels for each
#'   cell. Row meta-data can be accessed with \code{\link[SummarizedExperiment]{rowData}}.
#' 
#' 
#' @importFrom FlowSOM ReadInput BuildSOM BuildMST
#' @importFrom flowCore flowFrame
#' @importFrom SummarizedExperiment assays rowData colData 'rowData<-'
#' @importFrom grDevices pdf dev.off
#' 
#' @export
#' 
#' @seealso \code{\link{testDA}} \code{\link{testDE_med}} \code{\link{testDE_FDA}}
#'
#' @examples
#' # See full examples in testing functions: testDA, testDE_med, testDE_FDA
#' 
generateClusters <- function(d_se, 
                             cols_to_use = NULL, 
                             xdim = 20, ydim = 20, 
                             seed = NULL, 
                             plot = TRUE, path = ".", filename = "plot_MST.pdf", ...) {
  
  if (is.null(cols_to_use)) cols_to_use <- colData(d_se)$is_lineage
  
  if (!is.null(seed)) set.seed(seed)
  
  # FlowSOM requires input data as 'flowFrame' or 'flowSet'
  d_ff <- flowFrame(assays(d_se)[[1]])
  
  # note: not using FlowSOM 'meta-clustering' step
  runtime <- system.time({
    fsom <- ReadInput(d_ff, transform = FALSE, scale = FALSE); 
    fsom <- BuildSOM(fsom, colsToUse = cols_to_use, xdim = xdim, ydim = ydim); 
    fsom <- BuildMST(fsom)
  })
  
  message("FlowSOM clustering completed in ", round(runtime[["elapsed"]], 1), " seconds")
  
  # [to do: move to separate plotting function; generalize arguments]
  # save minimum spanning tree (MST) plot
  if (plot) {
    pdf(file.path(path, filename), width = 7, height = 7)
    plot(fsom$MST$l, cex = (fsom$MST$size - 7.5) / 3, pch = 19)
    dev.off()
  }
  
  # cluster labels
  clus <- fsom$map$mapping[, 1]
  
  # convert to factor (levels in ascending order)
  clus <- factor(clus, levels = sort(unique(clus)))
  
  # add cluster labels to row meta-data
  rowData(d_se)$cluster <- clus
  
  d_se
}


