#' Run 'diffcyt' pipeline
#' 
#' Wrapper function to run complete 'diffcyt' pipeline
#' 
#' This wrapper function runs the complete 'diffcyt' analysis pipeline, by calling each
#' individual function in the correct sequence.
#' 
#' For more details about the functions for each step in the pipeline, see the package
#' vignette and the individual function help pages. Running the individual functions may
#' provide additional flexibility, especially for complex analyses.
#' 
#' Minimum required arguments are:
#' 
#' \itemize{
#' \item \code{d_input}
#' \item \code{experiment_info}
#' \item \code{marker_info}
#' \item either \code{design} or \code{formula} (depending on the differential testing
#' method used)
#' \item \code{contrast}
#' \item \code{analysis_type}
#' }
#' 
#' 
#' @param d_input Input data. Must be a \linkS4class{flowSet} or list of
#'   \code{flowFrames}, \code{data.frames}, or matrices as input (one \code{flowFrame} or
#'   list item per sample). See \code{\link{prepareData}}.
#' 
#' @param experiment_info Data frame of experiment information, for example sample IDs and
#'   group IDs. Must contain a column named \code{sample_id}. See
#'   \code{\link{prepareData}}.
#' 
#' @param marker_info Data frame of marker information for each column of data. This
#'   should contain columns named \code{marker_name} and \code{marker_class}. The columns
#'   contain: (i) marker names (and any other column names); and (ii) a factor indicating
#'   the marker class for each column (with entries \code{"cell_type"},
#'   \code{"cell_state"}, or \code{"none"}). See \code{\link{prepareData}}.
#' 
#' @param design Design matrix, created with \code{\link{createDesignMatrix}}. See
#'   \code{\link{createDesignMatrix}}.
#' 
#' @param formula Model formula object, created with \code{\link{createFormula}}. See
#'   \code{\link{createFormula}}.
#' 
#' @param contrast Contrast matrix, created with \code{\link{createContrast}}. See
#'   \code{\link{createContrast}}.
#' 
#' @param analysis_type Type of differential analysis to perform: differential abundance
#'   (DA) of cell populations, or differential states (DS) within cell populations.
#'   Options are \code{"DA"} and \code{"DS"}. See \code{\link{testDA_edgeR}},
#'   \code{\link{testDA_voom}}, \code{\link{testDA_GLMM}}, \code{\link{testDS_limma}}, or
#'   \code{\link{testDS_LMM}}.
#' 
#' @param method_DA Method to use for calculating differential abundance (DA) tests.
#'   Options are \code{"diffcyt-DA-edgeR"}, \code{"diffcyt-DA-voom"}, and
#'   \code{"diffcyt-DA-GLMM"}. Default = \code{"diffcyt-DA-edgeR"}. See
#'   \code{\link{testDA_edgeR}}, \code{\link{testDA_voom}}, or \code{\link{testDA_GLMM}}.
#' 
#' @param method_DS Method to use for calculating differential state (DS) tests. Options
#'   are \code{"diffcyt-DS-limma"} and \code{"diffcyt-DS-LMM"}. Default =
#'   \code{"diffcyt-DS-limma"}. See \code{\link{testDS_limma}} or
#'   \code{\link{testDS_LMM}}.
#' 
#' @param cols_to_include Logical vector indicating which columns to include from the
#'   input data. Default = all columns. See \code{\link{prepareData}}.
#' 
#' @param subsampling Whether to use random subsampling to select an equal number of cells
#'   from each sample. Default = FALSE. See \code{\link{prepareData}}.
#' 
#' @param n_sub Number of cells to select from each sample by random subsampling, if
#'   \code{subsampling = TRUE}. Default = number of cells in smallest sample. See
#'   \code{\link{prepareData}}.
#' 
#' @param seed_sub Random seed for subsampling. Set to an integer value to generate
#'   reproducible results. Default = \code{NULL}. See \code{\link{prepareData}}.
#' 
#' @param cofactor Cofactor parameter for 'arcsinh' transform. Default = 5, which is
#'   appropriate for mass cytometry (CyTOF) data. For fluorescence flow cytometry, we
#'   recommend cofactor = 150 instead. See \code{\link{transformData}}.
#' 
#' @param cols_clustering Columns to use for clustering. Default = \code{NULL}, in which
#'   case markers identified as 'cell type' markers (with entries \code{"cell_type"}) in
#'   the vector \code{marker_class} in the column meta-data of \code{d_se} will be used.
#'   See \code{\link{generateClusters}}.
#' 
#' @param xdim Horizontal length of grid for self-organizing map for FlowSOM clustering
#'   (number of clusters = \code{xdim} * \code{ydim}). Default = 10 (i.e. 100 clusters).
#'   See \code{\link{generateClusters}}.
#' 
#' @param ydim Vertical length of grid for self-organizing map for FlowSOM clustering
#'   (number of clusters = \code{xdim} * \code{ydim}). Default = 10 (i.e. 100 clusters).
#'   See \code{\link{generateClusters}}.
#' 
#' @param meta_clustering Whether to include FlowSOM 'meta-clustering' step. Default =
#'   \code{FALSE}. See \code{\link{generateClusters}}.
#' 
#' @param meta_k Number of meta-clusters for FlowSOM, if \code{meta-clustering = TRUE}.
#'   Default = 40. See \code{\link{generateClusters}}.
#' 
#' @param seed_clustering Random seed for clustering. Set to an integer value to generate
#'   reproducible results. Default = \code{NULL}. See \code{\link{generateClusters}}.
#' 
#' @param min_cells Filtering parameter. Default = 3. Clusters are kept for differential
#'   testing if they have at least \code{min_cells} cells in at least \code{min_samples}
#'   samples. See \code{\link{testDA_edgeR}}, \code{\link{testDA_voom}},
#'   \code{\link{testDA_GLMM}}, \code{\link{testDS_limma}}, or \code{\link{testDS_LMM}}.
#' 
#' @param min_samples Filtering parameter. Default = \code{number of samples / 2}, which
#'   is appropriate for two-group comparisons (of equal size). Clusters are kept for
#'   differential testing if they have at least \code{min_cells} cells in at least
#'   \code{min_samples} samples. See \code{\link{testDA_edgeR}},
#'   \code{\link{testDA_voom}}, \code{\link{testDA_GLMM}}, \code{\link{testDS_limma}}, or
#'   \code{\link{testDS_LMM}}.
#' 
#' @param normalize Whether to include optional normalization factors to adjust for
#'   composition effects. Default = FALSE. See \code{\link{testDA_edgeR}},
#'   \code{\link{testDA_voom}}, or \code{\link{testDA_GLMM}}.
#' 
#' @param norm_factors Normalization factors to use, if \code{normalize = TRUE}. Default =
#'   \code{"TMM"}, in which case normalization factors are calculated automatically using
#'   the 'trimmed mean of M-values' (TMM) method from the \code{edgeR} package.
#'   Alternatively, a vector of values can be provided. (Note that other normalization
#'   methods from \code{edgeR} are not used.) See \code{\link{testDA_edgeR}},
#'   \code{\link{testDA_voom}}, or \code{\link{testDA_GLMM}}.
#' 
#' @param block_id (Optional) Vector or factor of block IDs (e.g. patient IDs) for paired
#'   experimental designs, to be included as random effects (for method \code{testDA_voom}
#'   or \code{testDS_limma}). If provided, the block IDs will be included as random
#'   effects using the \code{limma} \code{duplicateCorrelation} methodology.
#'   Alternatively, block IDs can be included as fixed effects in the design matrix
#'   (\code{\link{createDesignMatrix}}). See \code{\link{testDA_voom}} or
#'   \code{\link{testDS_limma}}.
#' 
#' @param plot Whether to save diagnostic plots (for method \code{testDA_voom} or
#'   \code{testDS_limma}). Default = TRUE. See \code{\link{testDA_voom}} or
#'   \code{\link{testDS_limma}}.
#' 
#' @param path Path for diagnostic plots (for method \code{testDA_voom} or
#'   \code{testDS_limma}). Default = current working directory. See
#'   \code{\link{testDA_voom}} or \code{\link{testDS_limma}}.
#' 
#' @param markers_to_test (Optional) Logical vector specifying which markers to test for
#'   differential expression (from the set of markers stored in the \code{assays} of
#'   \code{d_medians}; for method \code{\link{testDS_limma}} or \code{\link{testDS_LMM}}).
#'   Default = all 'cell state' markers, which are identified by the logical vector
#'   \code{id_state_markers} stored in the meta-data of \code{d_medians}. See
#'   \code{\link{testDS_limma}} or \code{\link{testDS_LMM}}.
#' 
#' 
#' @return Returns a list containing the results object \code{res}, as well as the data
#'   objects \code{d_se}, \code{d_counts}, \code{d_medians},
#'   \code{d_medians_by_cluster_marker}, and \code{d_medians_by_sample_marker}. The
#'   structure of \code{res} depends on the differential testing method used. See
#'   \code{\link{testDA_edgeR}}, \code{\link{testDA_voom}}, \code{\link{testDA_GLMM}},
#'   \code{\link{testDS_limma}}, or \code{\link{testDS_LMM}}.
#' 
#' 
#' @importFrom SummarizedExperiment colData
#' 
#' @export
#' 
#' @examples
#' # See the package vignette for a full workflow example showing both types of 
#' # differential analysis (DA and DS), and demonstrating each function in the 'diffcyt' 
#' # pipeline.
#' 
#' # Function to create random data (one sample)
#' d_random <- function(n = 20000, mean = 0, sd = 1, ncol = 20, cofactor = 5) {
#'   sinh(matrix(rnorm(n, mean, sd), ncol = ncol)) * cofactor
#' }
#' # Create random data (without differential signal)
#' set.seed(123)
#' d_input <- list(
#'   sample1 = d_random(), 
#'   sample2 = d_random(), 
#'   sample3 = d_random(), 
#'   sample4 = d_random()
#' )
#' # Add differential signal (for some cells and markers in one group)
#' ix_rows <- 901:1000
#' ix_cols <- c(6:10, 16:20)
#' d_input[[3]][ix_rows, ix_cols] <- d_random(n = 1000, mean = 3, ncol = 10)
#' d_input[[4]][ix_rows, ix_cols] <- d_random(n = 1000, mean = 3, ncol = 10)
#' 
#' experiment_info <- data.frame(
#'   sample_id = factor(paste0("sample", 1:4)), 
#'   group_id = factor(c("group1", "group1", "group2", "group2")), 
#'   stringsAsFactors = FALSE
#' )
#' 
#' marker_info <- data.frame(
#'   marker_name = paste0("marker", sprintf("%02d", 1:20)), 
#'   marker_class = factor(c(rep("cell_type", 10), rep("cell_state", 10)), 
#'                         levels = c("cell_type", "cell_state", "none")), 
#'   stringsAsFactors = FALSE
#' )
#' 
#' # Create design matrix
#' design <- createDesignMatrix(experiment_info, cols_design = 2)
#' # Create contrast matrix
#' contrast <- createContrast(c(0, 1))
#' 
#' # Test for differential abundance (DA) of clusters (using default method 'diffcyt-DA-edgeR')
#' out_DA <- diffcyt(d_input, experiment_info, marker_info, 
#'                   design = design, contrast = contrast, 
#'                   analysis_type = "DA", method_DA = "diffcyt-DA-edgeR")
#' 
#' # Test for differential states (DS) within clusters (using default method 'diffcyt-DS-limma')
#' out_DS <- diffcyt(d_input, experiment_info, marker_info, 
#'                   design = design, contrast = contrast, 
#'                   analysis_type = "DS", method_DS = "diffcyt-DS-limma", 
#'                   plot = FALSE)
#' 
#' # Display results for top DA clusters
#' topClusters(out_DA$res)
#' 
#' # Display results for top DS cluster-marker combinations
#' topClusters(out_DS$res)
#' 
#' # Plot heatmap for DA tests
#' plotHeatmap(out_DA, analysis_type = "DA")
#' 
#' # Plot heatmap for DS tests
#' plotHeatmap(out_DS, analysis_type = "DS")
#' 
diffcyt <- function(d_input, experiment_info, marker_info, design = NULL, formula = NULL, contrast, 
                    analysis_type = c("DA", "DS"), 
                    method_DA = c("diffcyt-DA-edgeR", "diffcyt-DA-voom", "diffcyt-DA-GLMM"), 
                    method_DS = c("diffcyt-DS-limma", "diffcyt-DS-LMM"), 
                    cols_to_include = NULL, 
                    subsampling = FALSE, n_sub = NULL, seed_sub = NULL, 
                    cofactor = 5, 
                    cols_clustering = NULL, xdim = 10, ydim = 10, 
                    meta_clustering = FALSE, meta_k = 40, seed_clustering = NULL, 
                    min_cells = 3, min_samples = NULL, 
                    normalize = FALSE, norm_factors = "TMM", 
                    block_id = NULL, 
                    plot = TRUE, path = ".", 
                    markers_to_test = NULL) {
  
  # arguments
  analysis_type <- match.arg(analysis_type)
  method_DA <- match.arg(method_DA)
  method_DS <- match.arg(method_DS)
  
  # prepare data, transformation, clustering
  d_se <- prepareData(d_input, experiment_info, marker_info, cols_to_include, subsampling, n_sub, seed_sub)
  d_se <- transformData(d_se, cofactor)
  d_se <- generateClusters(d_se, cols_clustering, xdim, ydim, meta_clustering, meta_k, seed_clustering)
  
  # calculate features
  d_counts <- calcCounts(d_se)
  d_medians <- calcMedians(d_se)
  
  # calculate additional features for plotting
  d_medians_by_cluster_marker <- calcMediansByClusterMarker(d_se)
  d_medians_by_sample_marker <- calcMediansBySampleMarker(d_se)
  
  # DA tests
  if (analysis_type == "DA" & method_DA == "diffcyt-DA-edgeR") {
    res <- testDA_edgeR(d_counts, design, contrast, min_cells, min_samples, normalize, norm_factors)
  }
  if (analysis_type == "DA" & method_DA == "diffcyt-DA-voom") {
    res <- testDA_voom(d_counts, design, contrast, block_id, min_cells, min_samples, normalize, norm_factors, plot, path)
  }
  if (analysis_type == "DA" & method_DA == "diffcyt-DA-GLMM") {
    res <- testDA_GLMM(d_counts, formula, contrast, min_cells, min_samples, normalize, norm_factors)
  }
  
  # DS tests
  if (analysis_type == "DS" & method_DS == "diffcyt-DS-limma") {
    res <- testDS_limma(d_counts, d_medians, design, contrast, block_id, markers_to_test, min_cells, min_samples, plot, path)
  }
  if (analysis_type == "DS" & method_DS == "diffcyt-DS-LMM") {
    res <- testDS_LMM(d_counts, d_medians, formula, contrast, markers_to_test, min_cells, min_samples)
  }
  
  # return results and data objects
  list(
    res = res, 
    d_se = d_se, 
    d_counts = d_counts, 
    d_medians = d_medians, 
    d_medians_by_cluster_marker = d_medians_by_cluster_marker, 
    d_medians_by_sample_marker = d_medians_by_sample_marker
  )
}


