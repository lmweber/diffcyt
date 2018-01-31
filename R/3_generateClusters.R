#' Generate clusters
#' 
#' Generate high-resolution clusters for \code{diffcyt} analysis
#' 
#' Performs clustering to group cells into clusters representing cell populations or
#' subsets, which can then be further analyzed for differential abundance of cell
#' populations or differential states within cell populations. By default, we use
#' high-resolution clustering or over-clustering (i.e. we generate a large number of small
#' clusters), which helps ensure that rare populations are adequately separated from
#' larger ones.
#' 
#' Data is assumed to be in the form of a
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} object generated with
#' \code{\link{prepareData}} and transformed with \code{\link{transformData}}.
#' 
#' The input data object \code{d_se} is assumed to contain a vector
#' \code{is_celltype_marker} in the column meta-data. This vector specifies the columns
#' that contain protein markers used to define cell types; these markers will be used for
#' clustering. For example, in immunological data, this may be the lineage markers. The
#' choice of cell type markers is an important design choice for the user, and will depend
#' on the underlying experimental design and research questions. It may be made based on
#' prior biological knowledge or using data-driven methods. For an example of a
#' data-driven method of marker ranking and selection, see Nowicka et al. (2017),
#' \emph{F1000Research}.
#' 
#' We use the \code{\link[FlowSOM]{FlowSOM}} clustering algorithm (Van Gassen et al. 2015,
#' \emph{Cytometry Part A}, available from Bioconductor) to generate the clusters. We
#' previously showed that \code{FlowSOM} gives very good clustering performance for
#' high-dimensional cytometry data, for both major and rare cell populations, and is
#' extremely fast (Weber and Robinson, 2016, \emph{Cytometry Part A}).
#' 
#' By default, the clustering is run at high resolution to give a large number of small
#' clusters (i.e. over-clustering). This is done by running only the initial
#' 'self-organizing map' clustering step in the \code{FlowSOM} algorithm, i.e. without the
#' final 'meta-clustering' step. This ensures that small or rare populations are
#' adequately separated from larger populations, which is crucial for detecting
#' differential abundance or differential states for extremely rare populations.
#' 
#' The minimum spanning tree (MST) object from \code{\link[FlowSOM]{BuildMST}} is stored
#' in the experiment \code{metadata} slot in the
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} object \code{d_se}, and can be
#' accessed with \code{metadata(d_se)$MST}.
#' 
#' 
#' @param d_se Transformed input data, from \code{\link{prepareData}} and
#'   \code{\link{transformData}}.
#' 
#' @param cols_to_use Columns to use for clustering. Default = \code{NULL}, in which case
#'   the markers specified by \code{is_celltype_marker} in the column meta-data of
#'   \code{d_se} will be used.
#' 
#' @param xdim Horizontal length of grid for self-organizing map for FlowSOM clustering
#'   (number of clusters = \code{xdim} * \code{ydim}). Default = 20 (i.e. 400 clusters).
#' 
#' @param ydim Vertical length of grid for self-organizing map for FlowSOM clustering
#'   (number of clusters = \code{xdim} * \code{ydim}). Default = 20 (i.e. 400 clusters).
#' 
#' @param meta_clustering Whether to include FlowSOM 'meta-clustering' step. Default =
#'   \code{FALSE}.
#' 
#' @param meta_k Number of meta-clusters for FlowSOM, if \code{meta-clustering = TRUE}.
#'   Default = 40.
#' 
#' @param seed Random seed. Set to an integer value to generate reproducible results.
#'   Default = \code{NULL}.
#' 
#' @param ... Other parameters to pass to the FlowSOM clustering algorithm (through the
#'   function \code{\link[FlowSOM]{BuildSOM}}).
#' 
#' 
#' @return \code{d_se}: Returns the
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}} input object, with cluster
#'   labels for each cell stored in an additional column of row meta-data. Row meta-data
#'   can be accessed with \code{\link[SummarizedExperiment]{rowData}}. The minimum
#'   spanning tree (MST) object is also stored in the \code{metadata} slot, and can be
#'   accessed with \code{metadata(d_se)$MST}.
#' 
#' 
#' @importFrom FlowSOM ReadInput BuildSOM BuildMST metaClustering_consensus
#' @importFrom flowCore flowFrame
#' @importFrom SummarizedExperiment assays rowData colData 'rowData<-'
#' @importFrom S4Vectors metadata 'metadata<-'
#' @importFrom grDevices pdf dev.off
#' 
#' @export
#' 
#' @examples
#' # A full workflow example demonstrating the use of each function in the 'diffcyt'
#' # pipeline on an experimental data set is available in the package vignette.
#' 
generateClusters <- function(d_se, cols_to_use = NULL, 
                             xdim = 20, ydim = 20, 
                             meta_clustering = FALSE, meta_k = 40, 
                             seed = NULL, ...) {
  
  if (is.null(cols_to_use)) cols_to_use <- colData(d_se)$is_celltype_marker
  
  # note: FlowSOM requires input data as 'flowFrame' or 'flowSet'
  d_ff <- flowFrame(assays(d_se)[[1]])
  
  runtime <- system.time({
    # FlowSOM: pre meta-clustering
    if (!is.null(seed)) set.seed(seed); 
    fsom <- ReadInput(d_ff, transform = FALSE, scale = FALSE); 
    fsom <- suppressMessages(BuildSOM(fsom, colsToUse = cols_to_use, xdim = xdim, ydim = ydim, ...)); 
    fsom <- suppressMessages(BuildMST(fsom))
    
    if (meta_clustering) {
      # FlowSOM: meta-clustering
      # (note: seed for meta-clustering must be provided via argument)
      meta <- metaClustering_consensus(fsom$map$codes, k = meta_k, seed = seed)
    }
  })
  
  message("FlowSOM clustering completed in ", round(runtime[["elapsed"]], 1), " seconds")
  
  # extract cluster labels
  clus_pre <- fsom$map$mapping[, 1]
  if (meta_clustering) {
    clus <- meta[clus_pre]
  } else {
    clus <- clus_pre
  }
  
  # convert to factor (with levels in ascending order)
  if (meta_clustering) {
    n_clus <- meta_k
  } else {
    n_clus <- xdim * ydim
  }
  
  clus <- factor(clus, levels = seq_len(n_clus))  # includes levels for any missing/empty clusters
  
  # store cluster labels in row meta-data
  rowData(d_se)$cluster <- clus
  
  # store MST object in experiment 'metadata' slot
  metadata(d_se)$MST <- fsom$MST
  
  d_se
}


