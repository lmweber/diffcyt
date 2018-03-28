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
#' \item \code{sample_info}
#' \item \code{marker_info}
#' \item either \code{design} or \code{formula}
#' \item \code{contrast}
#' \item \code{analysis_type}
#' \item either \code{method_DA} or \code{method_DS}
#' }
#' 
#' 
#' @param d_input Input data. Must be a list or \code{\link[flowCore]{flowSet}} (one list
#'   item or \code{\link[flowCore]{flowFrame}} per sample). See \code{\link{prepareData}}.
#' 
#' @param sample_info Data frame of sample information, for example sample IDs and group
#'   IDs. Must contain a column named \code{sample_IDs}. See \code{\link{prepareData}}.
#' 
#' @param marker_info Data frame of marker information for each column. This should
#'   contain columns named \code{marker_names}, \code{is_marker}, \code{is_type_marker},
#'   and \code{is_state_marker}. The first column must contain marker names or column
#'   names; the remaining columns are logical vectors indicating whether each column in
#'   the input data is (i) a protein marker, (ii) a 'cell type' marker, and (iii) a 'cell
#'   state' marker. See \code{\link{prepareData}}.
#' 
#' @param design Design matrix, created with \code{\link{createDesignMatrix}}. See
#'   \code{\link{createDesignMatrix}}.
#' 
#' @param formula Model formula object, created with \code{\link{createFormula}}. This
#'   should be a list containing three elements: \code{formula}, \code{data}, and
#'   \code{random_terms}: the model formula, data frame of corresponding variables, and
#'   variable indicating whether the model formula contains any random effect terms. See
#'   \code{\link{createFormula}}.
#' 
#' @param contrast Contrast matrix, created with \code{\link{createContrast}}. See
#'   \code{\link{createContrast}}.
#' 
#' @param analysis_type Type of differential discovery analysis to perform: differential
#'   abundance (DA) of cell populations; or differential states (DS) within cell
#'   populations. Options are \code{"DA"} and \code{"DS"}. See \code{\link{testDA_edgeR}},
#'   \code{\link{testDA_voom}}, \code{\link{testDA_GLMM}}, \code{\link{testDS_limma}}, or
#'   \code{\link{testDS_LMM}}.
#' 
#' @param method_DA Differential testing method to use, for DA tests. Options are
#'   \code{"diffcyt-DA-edgeR"}, \code{"diffcyt-DA-voom"}, and \code{"diffcyt-DA-GLMM"}.
#'   Default = \code{"diffcyt-DA-edgeR"}. See \code{\link{testDA_edgeR}},
#'   \code{\link{testDA_voom}}, or \code{\link{testDA_GLMM}}.
#' 
#' @param method_DS Differential testing method to use, for DS tests. Options are
#'   \code{"diffcyt-DS-limma"} and \code{"diffcyt-DS-LMM"}. Default =
#'   \code{"diffcyt-DS-limma"}. See \code{\link{testDS_limma}} or
#'   \code{\link{testDS_LMM}}.
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
#' @param cols_to_use Columns to use for clustering. Default = \code{NULL}, in which case
#'   the markers specified by \code{is_type_marker} in the column meta-data of \code{d_se}
#'   will be used. See \code{\link{generateClusters}}.
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
#'   composition effects (see details). Default = FALSE. See \code{\link{testDA_edgeR}},
#'   \code{\link{testDA_voom}}, or \code{\link{testDA_GLMM}}.
#' 
#' @param norm_factors Normalization factors to use, if \code{normalize = TRUE}. Default =
#'   \code{"TMM"}, in which case normalization factors are calculated automatically using
#'   the 'trimmed mean of M-values' (TMM) method from the \code{edgeR} package. See
#'   \code{\link{testDA_edgeR}}, \code{\link{testDA_voom}}, or \code{\link{testDA_GLMM}}.
#' 
#' @param block_IDs (Optional) Vector or factor of block IDs (e.g. patient IDs) for paired
#'   experimental designs, to be included as random effects. If provided, the block IDs
#'   will be included as random effects using the \code{limma}
#'   \code{\link[limma]{duplicateCorrelation}} methodology. Alternatively, block IDs can
#'   be included as fixed effects in the design matrix (\code{\link{createDesignMatrix}}).
#'   See details. See \code{\link{testDA_voom}} or \code{\link{testDS_limma}}.
#' 
#' @param plot Whether to save diagnostic plots. Default = TRUE. See
#'   \code{\link{testDA_voom}} or \code{\link{testDS_limma}}.
#' 
#' @param path Path for diagnostic plots. Default = current working directory. See
#'   \code{\link{testDA_voom}} or \code{\link{testDS_limma}}.
#' 
#' 
#' @return Returns the results object from the differential testing method. The structure
#'   of the object depends on which method was used. See \code{\link{testDA_edgeR}},
#'   \code{\link{testDA_voom}}, \code{\link{testDA_GLMM}}, \code{\link{testDS_limma}}, or
#'   \code{\link{testDS_LMM}}.
#' 
#' 
#' @importFrom SummarizedExperiment colData
#' 
#' @export
#' 
#' @examples
#' # See the package vignette for a full workflow example demonstrating each type of
#' # differential discovery analysis (DA and DS), and explaining each function in the
#' # 'diffcyt' pipeline.
#' 
#' # Create some random data (without differential signal)
#' cofactor <- 5
#' set.seed(123)
#' d_input <- list(
#'   sample1 = sinh(matrix(rnorm(20000, mean = 0, sd = 1), ncol = 20)) * cofactor, 
#'   sample2 = sinh(matrix(rnorm(20000, mean = 0, sd = 1), ncol = 20)) * cofactor, 
#'   sample3 = sinh(matrix(rnorm(20000, mean = 0, sd = 1), ncol = 20)) * cofactor, 
#'   sample4 = sinh(matrix(rnorm(20000, mean = 0, sd = 1), ncol = 20)) * cofactor
#' )
#' # Add differential signal (for some cells and markers in one group)
#' ix_rows <- 901:1000
#' ix_cols <- c(6:10, 16:20)
#' d_input[[3]][ix_rows, ix_cols] <- sinh(matrix(rnorm(1000, mean = 2, sd = 1), ncol = 10)) * cofactor
#' d_input[[4]][ix_rows, ix_cols] <- sinh(matrix(rnorm(1000, mean = 2, sd = 1), ncol = 10)) * cofactor
#' 
#' sample_info <- data.frame(
#'   sample_IDs = paste0("sample", 1:4), 
#'   group_IDs = factor(c("group1", "group1", "group2", "group2"))
#' )
#' 
#' marker_info <- data.frame(
#'   marker_names = paste0("marker", 1:20), 
#'   is_marker = rep(TRUE, 20), 
#'   is_type_marker = c(rep(TRUE, 10), rep(FALSE, 10)), 
#'   is_state_marker = c(rep(FALSE, 10), rep(TRUE, 10))
#' )
#' 
#' # Create design matrix
#' design <- createDesignMatrix(sample_info, cols_include = 2)
#' # Create contrast matrix
#' contrast <- createContrast(c(0, 1))
#' 
#' # Test for differential abundance (DA) of clusters (using default method 'diffcyt-DA-edgeR')
#' res_DA <- diffcyt(d_input, sample_info, marker_info, 
#'                   design = design, contrast = contrast, 
#'                   analysis_type = "DA", method_DA = "diffcyt-DA-edgeR")
#' 
#' # Test for differential states (DS) within clusters (using default method 'diffcyt-DS-limma')
#' res_DS <- diffcyt(d_input, sample_info, marker_info, 
#'                   design = design, contrast = contrast, 
#'                   analysis_type = "DS", method_DS = "diffcyt-DS-limma", 
#'                   plot = FALSE)
#' 
diffcyt <- function(d_input, sample_info, marker_info, design = NULL, formula = NULL, contrast, 
                    analysis_type = c("DA", "DS"), 
                    method_DA = c("diffcyt-DA-edgeR", "diffcyt-DA-voom", "diffcyt-DA-GLMM"), 
                    method_DS = c("diffcyt-DS-limma", "diffcyt-DS-LMM"), 
                    subsampling = FALSE, n_sub = NULL, seed_sub = NULL, 
                    cofactor = 5, 
                    cols_to_use = NULL, xdim = 10, ydim = 10, 
                    meta_clustering = FALSE, meta_k = 40, seed_clustering = NULL, 
                    min_cells = 3, min_samples = NULL, normalize = FALSE, norm_factors = "TMM", 
                    block_IDs = NULL, plot = TRUE, path = ".") {
  
  # arguments
  analysis_type <- match.arg(analysis_type)
  method_DA <- match.arg(method_DA)
  method_DS <- match.arg(method_DS)
  
  # prepare data, transformation, clustering
  d_se <- prepareData(d_input, sample_info, marker_info, subsampling, n_sub, seed_sub)
  d_se <- transformData(d_se, cofactor)
  d_se <- generateClusters(d_se, cols_to_use, xdim, ydim, meta_clustering, meta_k, seed_clustering)
  
  # calculate features
  d_counts <- calcCounts(d_se)
  d_medians <- calcMedians(d_se)
  
  # calculate tests
  if (analysis_type == "DA" & method_DA == "diffcyt-DA-edgeR") {
    res <- testDA_edgeR(d_counts, design, contrast, min_cells, min_samples, normalize, norm_factors)
  }
  if (analysis_type == "DA" & method_DA == "diffcyt-DA-voom") {
    res <- testDA_voom(d_counts, design, contrast, block_IDs, min_cells, min_samples, normalize, norm_factors, plot, path)
  }
  if (analysis_type == "DA" & method_DA == "diffcyt-DA-GLMM") {
    res <- testDA_GLMM(d_counts, formula, contrast, min_cells, min_samples, normalize, norm_factors)
  }
  if (analysis_type == "DS" & method_DS == "diffcyt-DS-limma") {
    res <- testDS_limma(d_counts, d_medians, design, contrast, block_IDs, min_cells, min_samples, plot, path)
  }
  if (analysis_type == "DS" & method_DS == "diffcyt-DS-LMM") {
    res <- testDS_LMM(d_counts, d_medians, formula, contrast, min_cells, min_samples)
  }
  
  # return test results
  res
}


