#' Generate clusters
#' 
#' Generate clusters for 'diffcyt' analysis
#' 
#' Runs a clustering algorithm to divide cells into a large number of small clusters, 
#' which will then be further analyzed for differential abundance or differential 
#' expression of functional markers.
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
#' In most cases, only lineage markers should be used for clustering.
#' 
#' Note that data should be transformed prior to clustering, which can be done with the
#' function \code{\link{transformData}}.
#' 
#' 
#' @param d_transf Transformed input data. Data should be transformed prior to clustering
#'   (see \code{\link{transformData}}). Input data must be in the form of a 
#'   \code{\link[flowCore]{flowSet}} object from the \code{\link[flowCore]{flowCore}} 
#'   package, which can be checked with \code{\link{checkInputData}}.
#' 
#' @param cols_to_use Columns to use for clustering. This will depend on the data set, 
#'   but in most cases will be the set of lineage markers. Default = NULL, in which case
#'   all columns will be used.
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
#' @param ... Other parameters to pass to FlowSOM (see \code{\link[FlowSOM]{BuildSOM}} 
#'   and \code{\link[FlowSOM]{SOM}} for details.)
#' 
#' 
#' @return clus List containing cluster labels for each cell (one list item per sample).
#' 
#' 
#' @importFrom FlowSOM ReadInput BuildSOM BuildMST
#' @importFrom flowCore pData
#' @importFrom grDevices pdf dev.off
#' 
#' @export
#' 
#' @examples
#' # need to create a small example data set for examples
generateClusters <- function(d_transf, 
                             cols_to_use = NULL, 
                             xdim = 20, ydim = 20, 
                             seed = NULL, plot = TRUE, ...) {
  
  if (is.null(cols_to_use)) cols_to_use <- 1:ncol(d_transf[[1]])
  
  if (!is.null(seed)) set.seed(seed)
  
  # note: not using meta-clustering
  runtime <- system.time({
    fsom <- FlowSOM::ReadInput(d_transf, transform = FALSE, scale = FALSE); 
    fsom <- FlowSOM::BuildSOM(fsom, colsToUse = cols_to_use, xdim = xdim, ydim = ydim); 
    fsom <- FlowSOM::BuildMST(fsom)
  })
  
  message("FlowSOM clustering completed in ", round(runtime[["elapsed"]], 1), " seconds")
  
  # [to do: move plotting to a separate function]
  if (plot) {
    pdf("plot_MST.pdf", width = 7, height = 7)
    plot(fsom$MST$l, cex = (fsom$MST$size - 7.5) / 3, pch = 19)
    dev.off()
  }
  
  # cluster labels
  clus_all <- fsom$map$mapping[,1]
  
  # split cluster labels by sample
  nm <- pData(d_transf)$name
  n_cells_smp <- sapply(as(d_transf, "list"), nrow)
  stopifnot(all(names(n_cells_smp) == nm))
  
  smp <- rep(nm, n_cells_smp)
  smp <- factor(smp, levels = nm)  # in case levels in non-alphabetical order
  
  clus <- split(clus_all, smp)
  
  clus
}


