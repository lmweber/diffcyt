#' Test for differential abundance: method 'diffcyt-DA-edgeR'
#' 
#' Calculate tests for differential abundance of cell populations using method
#' 'diffcyt-DA-edgeR'
#' 
#' Calculates tests for differential abundance of clusters, using functions from the
#' \code{edgeR} package.
#' 
#' This method uses the \code{edgeR} package (Robinson et al. 2010, \emph{Bioinformatics};
#' McCarthy et al. 2012, \emph{Nucleic Acids Research}) to fit models and calculate
#' moderated tests at the cluster level. Moderated tests improve statistical power by
#' sharing information on variability (i.e. variance across samples for a single cluster)
#' between clusters. We use the option \code{trend.method = "none"} to calculate
#' dispersions, since the dispersion-mean relationship typically does not resemble
#' RNA-sequencing data; see \code{edgeR} User's Guide. The statistical methods implemented
#' in the \code{edgeR} package were originally designed for the analysis of gene
#' expression data such as RNA-sequencing counts. Here, we apply these methods to cluster
#' cell counts.
#' 
#' The experimental design must be specified using a design matrix, which can be created
#' with \code{\link{createDesignMatrix}}. Flexible experimental designs are possible,
#' including blocking (e.g. paired designs), batch effects, and continuous covariates. See
#' \code{\link{createDesignMatrix}} for more details.
#' 
#' The contrast matrix specifying the contrast of interest can be created with
#' \code{\link{createContrast}}. See \code{\link{createContrast}} for more details.
#' 
#' Filtering: Clusters are kept for differential testing if they have at least
#' \code{min_cells} cells in at least \code{min_samples} samples. This removes clusters
#' with very low cell counts across conditions, to improve power.
#' 
#' Normalization for the total number of cells per sample (library sizes) and total number
#' of cells per cluster is automatically performed by the \code{edgeR} functions. Optional
#' normalization factors can also be included to adjust for composition effects in the
#' total number of counts per sample (library sizes). For example, if one cell population
#' is more abundant in one condition, while all other populations have equal numbers of
#' cells across conditions, then this effectively reduces the \emph{relative} abundance of
#' the unchanged populations in the first condition, creating false positive differential
#' abundance signals for these populations. Normalization factors can be provided directly
#' as a vector of values representing relative total abundances per sample (where values
#' >1 indicate extra cells in a sample; note this is the inverse of normalization factors
#' as defined in the \code{edgeR} package). Alternatively, normalization factors can be
#' calculated automatically using the 'trimmed mean of M-values' (TMM) method from the
#' \code{edgeR} package (Robinson and Oshlack, 2010). The TMM method assumes that most
#' populations are not differentially abundant; see the \code{edgeR} User's Guide for more
#' details.
#' 
#' 
#' @param d_counts \linkS4class{SummarizedExperiment} object containing cluster cell
#'   counts, from \code{\link{calcCounts}}.
#' 
#' @param design Design matrix, created with \code{\link{createDesignMatrix}}. See
#'   \code{\link{createDesignMatrix}} for details.
#' 
#' @param contrast Contrast matrix, created with \code{\link{createContrast}}. See
#'   \code{\link{createContrast}} for details.
#' 
#' @param min_cells Filtering parameter. Default = 3. Clusters are kept for differential
#'   testing if they have at least \code{min_cells} cells in at least \code{min_samples}
#'   samples.
#' 
#' @param min_samples Filtering parameter. Default = \code{number of samples / 2}, which
#'   is appropriate for two-group comparisons (of equal size). Clusters are kept for
#'   differential testing if they have at least \code{min_cells} cells in at least
#'   \code{min_samples} samples.
#' 
#' @param normalize Whether to include optional normalization factors to adjust for
#'   composition effects (see details). Default = FALSE.
#' 
#' @param norm_factors Normalization factors to use, if \code{normalize = TRUE}. Default =
#'   \code{"TMM"}, in which case normalization factors are calculated automatically using
#'   the 'trimmed mean of M-values' (TMM) method from the \code{edgeR} package.
#'   Alternatively, a vector of values can be provided. (Note that other normalization
#'   methods from \code{edgeR} are not used.)
#' 
#' 
#' @return Returns a new \linkS4class{SummarizedExperiment} object, with differential test
#'   results stored in the \code{rowData} slot. Results include raw p-values and adjusted
#'   p-values from the \code{edgeR} moderated tests, which can be used to rank clusters by
#'   evidence for differential abundance. The results can be accessed with the
#'   \code{rowData} accessor function.
#' 
#' 
#' @importFrom SummarizedExperiment assay rowData 'rowData<-' colData 'colData<-'
#' @importFrom edgeR calcNormFactors DGEList estimateDisp glmFit glmLRT topTags
#' @importFrom methods as is
#' 
#' @export
#' 
#' @examples
#' # For a full workflow example demonstrating the use of each function in the 'diffcyt'
#' # pipeline, see the package vignette.
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
#' # Prepare data
#' d_se <- prepareData(d_input, experiment_info, marker_info)
#' 
#' # Transform data
#' d_se <- transformData(d_se)
#' 
#' # Generate clusters
#' d_se <- generateClusters(d_se)
#' 
#' # Calculate counts
#' d_counts <- calcCounts(d_se)
#' 
#' # Create design matrix
#' design <- createDesignMatrix(experiment_info, cols_design = 2)
#' # Create contrast matrix
#' contrast <- createContrast(c(0, 1))
#' 
#' # Test for differential abundance (DA) of clusters
#' res_DA <- testDA_edgeR(d_counts, design, contrast)
#' 
testDA_edgeR <- function(d_counts, design, contrast, 
                         min_cells = 3, min_samples = NULL, 
                         normalize = FALSE, norm_factors = "TMM") {
  
  if (is.null(min_samples)) {
    min_samples <- ncol(d_counts) / 2
  }
  
  counts <- assay(d_counts)
  cluster_id <- rowData(d_counts)$cluster_id
  
  # filtering: keep clusters with at least 'min_cells' cells in at least 'min_samples' samples
  tf <- counts >= min_cells
  ix_keep <- apply(tf, 1, function(r) sum(r) >= min_samples)
  
  counts <- counts[ix_keep, , drop = FALSE]
  cluster_id <- cluster_id[ix_keep]
  
  # edgeR pipeline
  
  # normalization factors
  if (normalize & norm_factors == "TMM") {
    norm_factors <- calcNormFactors(counts, method = "TMM")
  } else if (normalize & norm_factors != "TMM") {
    # edgeR and limma define normalization factors as inverse and require product = 1
    norm_factors <- 1 / norm_factors
    # using geometric mean
    norm_factors <- norm_factors / (prod(norm_factors) ^ (1 / length(norm_factors)))
  }
  
  # note: when using DGEList object, column sums are automatically used for library sizes
  if (normalize) {
    y <- DGEList(counts, norm.factors = norm_factors)
  } else {
    y <- DGEList(counts)
  }
  
  # estimate dispersions
  # (note: using 'trend.method = "none"')
  y <- estimateDisp(y, design, trend.method = "none")
  
  # fit models
  fit <- glmFit(y, design)
  
  # likelihood ratio tests
  lrt <- glmLRT(fit, contrast = contrast)
  
  # results
  top <- topTags(lrt, n = Inf, adjust.method = "BH", sort.by = "none")
  
  if (!all(rownames(top) == cluster_id)) {
    stop("cluster labels do not match")
  }
  
  # return results in 'rowData' of new 'SummarizedExperiment' object
  
  # fill in any missing rows (filtered clusters) with NAs
  row_data <- as.data.frame(matrix(as.numeric(NA), nrow = nlevels(cluster_id), ncol = ncol(top)))
  colnames(row_data) <- colnames(top)
  cluster_id_nm <- as.numeric(cluster_id)
  row_data[cluster_id_nm, ] <- top$table
  
  row_data <- cbind(data.frame(cluster_id = as.numeric(levels(cluster_id)), stringsAsFactors = FALSE), 
                    row_data)
  
  res <- d_counts
  
  rowData(res) <- row_data
  
  res
}


