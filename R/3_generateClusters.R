#' Generate mini-clusters
#' 
#' Generate mini-clusters for \code{diffcyt} analysis
#' 
#' Performs clustering to divide cells into subsets (clusters or populations), which can 
#' then be further analyzed for differential abundance, or differential expression within 
#' clusters. By default, a large number of small clusters is generated ('mini-clusters', 
#' i.e. over-clustering), to ensure that rare populations are adequately separated from 
#' larger ones.
#' 
#' Data is assumed to be in the form of a 
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} object generated with 
#' \code{\link{prepareData}}, and already transformed with \code{\link{transformData}}.
#' 
#' The input data object \code{d_se} is assumed to contain a vector 
#' \code{is_clustering_col}, which specifies the columns (protein markers) to be used for 
#' clustering. For example, in immunological data, this may be the lineage markers. 
#' Alternatively, for a more sophisticated method of marker ranking and selection for 
#' clustering, see the CyTOF data analysis workflow by Nowicka et al. (2017), 
#' \emph{F1000Research}.
#' 
#' We use the \code{\link[FlowSOM]{FlowSOM}} clustering algorithm (Van Gassen et al. 2015,
#' \emph{Cytometry Part A}, available from Bioconductor) to generate the clusters. We 
#' previously showed that \code{FlowSOM} gives very good clustering performance for 
#' high-dimensional cytometry data, for both major and rare cell populations, and is 
#' extremely fast (Weber and Robinson, 2016, \emph{Cytometry Part A}).
#' 
#' By default, the clustering is run at high resolution to give a large number of small 
#' clusters ('mini-clusters' or over-clustering). This is done by running only the initial
#' 'self-organizing map' clustering step in the \code{FlowSOM} algorithm, i.e. without the
#' final 'meta-clustering' step. This ensures that small or rare populations are 
#' adequately separated from larger populations, which is crucial for detecting 
#' differential abundance or differential expression signals from these populations.
#' 
#' The minimum spanning tree (MST) object from \code{\link[FlowSOM]{BuildMST}} (required 
#' for plotting) is stored in the experiment \code{metadata} slot in the 
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} object \code{d_se}, and can be
#' accessed with \code{metadata(d_se)$MST}.
#' 
#' 
#' @param d_se Transformed input data, from \code{\link{prepareData}} and 
#'   \code{\link{transformData}}.
#' 
#' @param cols_to_use Columns to use for clustering. Default = \code{NULL}, in which case 
#'   the markers specified by \code{is_clustering_col} in the column meta-data of 
#'   \code{d_se} will be used.
#' 
#' @param xdim Width of grid for self-organizing map for FlowSOM clustering (number of 
#'   clusters = \code{xdim} * \code{ydim}). Default = 30 (i.e. 900 clusters).
#' 
#' @param ydim Height of grid for self-organizing map for FlowSOM clustering (number of 
#'   clusters = \code{xdim} * \code{ydim}). Default = 30 (i.e. 900 clusters).
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
#' @return Returns the \code{\link[SummarizedExperiment]{SummarizedExperiment}} input 
#'   object, with cluster labels for each cell stored in an additional column of row 
#'   meta-data. Row meta-data can be accessed with 
#'   \code{\link[SummarizedExperiment]{rowData}}. The minimum spanning tree (MST) object 
#'   is stored in the experiment \code{metadata} slot, and can be accessed with 
#'   \code{metadata(d_se)$MST}.
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
#' @seealso \code{\link{testDA}} \code{\link{testDE_med}} \code{\link{testDE_FDA}} 
#'   \code{\link{testDE_KS}} \code{\link{testDE_LM}}
#'
#' @examples
#' # See full examples in testing functions: testDA, testDE_med, testDE_FDA, testDE_KS,
#' # testDE_LM
#' 
generateClusters <- function(d_se, cols_to_use = NULL, 
                             xdim = 30, ydim = 30, 
                             meta_clustering = FALSE, meta_k = 40, 
                             seed = NULL, ...) {
  
  if (is.null(cols_to_use)) cols_to_use <- colData(d_se)$is_clustering_col
  
  # note: FlowSOM requires input data as 'flowFrame' or 'flowSet'
  d_ff <- flowFrame(assays(d_se)[[1]])
  
  runtime <- system.time({
    # FlowSOM: pre meta-clustering
    # (note: FlowSOM seed is not reproducible across operating systems)
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
  
  clus <- factor(clus, levels = 1:n_clus)  # includes levels for missing clusters
  
  # store cluster labels in row meta-data
  rowData(d_se)$cluster <- clus
  
  # store MST object in experiment 'metadata' slot
  metadata(d_se)$MST <- fsom$MST
  
  d_se
}


