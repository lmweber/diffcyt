#' Generate clusters
#' 
#' Generate clusters for \code{diffcyt} analysis
#' 
#' Runs a clustering algorithm to divide cells into subsets (clusters), which will then be
#' further analyzed for differential abundance (frequencies) of clusters, or differential 
#' expression of functional markers.
#' 
#' Data is assumed to be in the form of a 
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}}, and already transformed with 
#' \code{\link{transformData}}.
#' 
#' By default, clustering is performed using lineage markers only. The column meta-data of
#' the input \code{SummarizedExperiment} object is assumed to contain a vector of logical 
#' entries (\code{is_lineage}) indicating which columns represent lineage markers.
#' 
#' We use the \code{\link[FlowSOM]{FlowSOM}} clustering algorithm (Van Gassen et al. 2015,
#' \emph{Cytometry Part A}) to generate the clusters. We previously showed that 
#' \code{FlowSOM} gives very good clustering performance for high-dimensional cytometry 
#' data, for both major and rare cell populations, and is extremely fast (Weber and 
#' Robinson, 2016, \emph{Cytometry Part A}).
#' 
#' By default, the clustering is run at high resolution to give a large number of small 
#' clusters (this is done by running only the initial 'self-organizing map' clustering 
#' step in the \code{FlowSOM} algorithm; i.e. without the additional 'meta-clustering'
#' step). This is appropriate when testing for differential expression of functional
#' markers, since it ensures that subtle differential expression signals in small subsets
#' are not missed.
#' 
#' When testing only for differential abundance (frequencies) of clusters, it may be more 
#' useful to return larger, more easily interpretabel clusters. This can be done by 
#' setting \code{meta_clustering = TRUE} to include the final 'meta-clustering' step from 
#' the \code{FlowSOM} algorithm (and optionally, specify the number of meta-clusters 
#' \code{meta_k}).
#' 
#' The minimum spanning tree (MST) object from \code{\link[FlowSOM]{BuildMST}} (required 
#' for plotting) is stored in the experiment \code{metadata} list in the \code{d_se}
#' object, and can be accessed with \code{metadata(d_se)$MST}.
#' 
#' 
#' @param d_se Transformed input data, from \code{\link{transformData}}.
#' 
#' @param cols_to_use Columns to use for clustering. Default = \code{NULL}, in which case 
#'   the lineage markers from \code{is_lineage} in the column meta-data of \code{d_se} 
#'   will be used.
#' 
#' @param xdim Width of grid for self-organizing map for FlowSOM clustering (number of 
#'   clusters = \code{xdim} * \code{ydim}). Default = 20 (i.e. 400 clusters).
#' 
#' @param ydim Height of grid for self-organizing map for FlowSOM clustering (number of 
#'   clusters = \code{xdim} * \code{ydim}). Default = 20 (i.e. 400 clusters).
#' 
#' @param meta_clustering Whether to include FlowSOM 'meta-clustering' step. Default = 
#'   \code{FALSE}.
#' 
#' @param meta_k Number of final meta-clusters for FlowSOM, if \code{meta-clustering = 
#'   TRUE}. Default = 40.
#' 
#' @param seed Random seed (set to an integer value for reproducible results). Default = 
#'   \code{NULL}.
#' 
#' @param ... Other parameters to pass to the FlowSOM clustering algorithm (through the 
#'   function \code{\link[FlowSOM]{BuildSOM}}).
#' 
#' 
#' @return Returns the \code{\link[SummarizedExperiment]{SummarizedExperiment}} input 
#'   object with an additional column of row meta-data, containing cluster labels for each
#'   cell. Row meta-data can be accessed with \code{\link[SummarizedExperiment]{rowData}}.
#'   The minimum spanning tree (MST) object is stored in the experiment \code{metadata}
#'   list, and can be accessed with \code{metadata(d_se)$MST}.
#' 
#' 
#' @importFrom FlowSOM ReadInput BuildSOM BuildMST
#' @importFrom ConsensusClusterPlus ConsensusClusterPlus
#' @importFrom flowCore flowFrame
#' @importFrom SummarizedExperiment assays rowData colData 'rowData<-'
#' @importFrom S4Vectors 'metadata<-'
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
                             xdim = 20, ydim = 20, meta_clustering = FALSE, meta_k = 40, 
                             seed = NULL, ...) {
  
  if (is.null(cols_to_use)) cols_to_use <- colData(d_se)$is_lineage
  
  if (!is.null(seed)) set.seed(seed)
  
  # note: FlowSOM requires input data as 'flowFrame' or 'flowSet'
  d_ff <- flowFrame(assays(d_se)[[1]])
  
  runtime <- system.time({
    # FlowSOM: pre meta-clustering
    fsom <- ReadInput(d_ff, transform = FALSE, scale = FALSE); 
    fsom <- BuildSOM(fsom, colsToUse = cols_to_use, xdim = xdim, ydim = ydim, ...); 
    fsom <- BuildMST(fsom)
    
    if (meta_clustering) {
      # FlowSOM: optional meta-clustering step (note: using ConsensusClusterPlus directly 
      # due to bug in random seed in 'FlowSOM::metaClustering_consensus()')
      meta <- ConsensusClusterPlus(t(fsom$map$codes), maxK = meta_k, seed = seed)
      meta <- meta[[meta_k]]$consensusClass
    }
  })
  
  message("FlowSOM clustering completed in ", round(runtime[["elapsed"]], 1), " seconds")
  
  # cluster labels
  clus_pre <- fsom$map$mapping[, 1]  # pre meta-clustering labels
  if (meta_clustering) {
    clus <- meta[clus_pre]  # meta-clustering labels
  } else {
    clus <- clus_pre
  }
  
  # convert to factor (levels in ascending order)
  clus <- factor(clus, levels = sort(unique(clus)))
  
  # add cluster labels to row meta-data
  rowData(d_se)$cluster <- clus
  
  # store MST object in experiment 'metadata' slot
  metadata(d_se)$MST <- fsom$MST
  
  d_se
}


