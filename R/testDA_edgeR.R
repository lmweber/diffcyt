#' Test for differential abundance: method 'diffcyt-DA-edgeR'
#' 
#' Calculate tests for differential abundance of cell populations using method
#' 'diffcyt-DA-edgeR'
#' 
#' Calculates tests for differential abundance of clusters, using functions from the
#' \code{\link{edgeR}} package.
#' 
#' This method uses the \code{\link{edgeR}} package (Robinson et al. 2010,
#' \emph{Bioinformatics}; McCarthy et al. 2012, \emph{Nucleic Acids Research}) to fit
#' models and calculate moderated tests at the cluster level. Moderated tests improve
#' statistical power by sharing information on variability (i.e. variance across samples
#' for a single cluster) between clusters. By default, we use the option
#' \code{trend.method = "none"} to calculate dispersions, since the dispersion-mean
#' relationship typically does not resemble RNA-sequencing data; see \code{edgeR} User's
#' Guide. The statistical methods implemented in the \code{edgeR} package were originally
#' designed for the analysis of gene expression data such as RNA-sequencing counts. Here,
#' we apply these methods to cluster cell counts.
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
#' cluster cell counts per sample. For example, in an extreme case, if several additional
#' clusters are present in only one condition, while all other clusters are approximately
#' equally abundant between conditions, then simply normalizing by the total number of
#' cells per sample will create a false positive differential abundance signal for the
#' non-differential clusters. (For a detailed explanation in the context of RNA sequencing
#' gene expression, see Robinson and Oshlack, 2010.) Normalization factors can be
#' calculated automatically using the 'trimmed mean of M-values' (TMM) method (Robinson
#' and Oshlack, 2010), implemented in the \code{edgeR} package (see also the \code{edgeR}
#' User's Guide for details). Alternatively, a vector of values can be provided (the
#' values should multiply to 1).
#' 
#' 
#' @param d_counts \code{\link{SummarizedExperiment}} object containing cluster cell
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
#'   Alternatively, a vector of values can be provided (the values should multiply to 1).
#' 
#' @param trend_method Method for estimating dispersion trend; passed to function
#'   \code{\link{estimateDisp}} from \code{edgeR} package. Default = "none". (See
#'   \code{\link{estimateDisp}} help file from \code{edgeR} package for other options.)
#' 
#' 
#' @return Returns a new \code{\link{SummarizedExperiment}} object, with differential test
#'   results stored in the \code{rowData} slot. Results include raw p-values
#'   (\code{p_val}) and adjusted p-values (\code{p_adj}) from the \code{edgeR} moderated
#'   tests, which can be used to rank clusters by evidence for differential abundance.
#'   Additional output columns from the \code{edgeR} tests are also included. The results
#'   can be accessed with the \code{\link{rowData}} accessor function.
#' 
#' 
#' @importFrom SummarizedExperiment assays rowData 'rowData<-' colData 'colData<-'
#' @importFrom edgeR calcNormFactors DGEList estimateDisp glmFit glmLRT topTags
#' @importFrom methods as is
#' 
#' @export
#' 
#' @examples
#' # For a complete workflow example demonstrating each step in the 'diffcyt' pipeline, 
#' # see the package vignette.
#' 
#' # Function to create random data (one sample)
#' d_random <- function(n = 20000, mean = 0, sd = 1, ncol = 20, cofactor = 5) {
#'   d <- sinh(matrix(rnorm(n, mean, sd), ncol = ncol)) * cofactor
#'   colnames(d) <- paste0("marker", sprintf("%02d", 1:ncol))
#'   d
#' }
#' 
#' # Create random data (without differential signal)
#' set.seed(123)
#' d_input <- list(
#'   sample1 = d_random(), 
#'   sample2 = d_random(), 
#'   sample3 = d_random(), 
#'   sample4 = d_random()
#' )
#' 
#' # Add differential abundance (DA) signal
#' ix_DA <- 801:900
#' ix_cols_type <- 1:10
#' d_input[[3]][ix_DA, ix_cols_type] <- d_random(n = 1000, mean = 2, ncol = 10)
#' d_input[[4]][ix_DA, ix_cols_type] <- d_random(n = 1000, mean = 2, ncol = 10)
#' 
#' experiment_info <- data.frame(
#'   sample_id = factor(paste0("sample", 1:4)), 
#'   group_id = factor(c("group1", "group1", "group2", "group2")), 
#'   stringsAsFactors = FALSE
#' )
#' 
#' marker_info <- data.frame(
#'   channel_name = paste0("channel", sprintf("%03d", 1:20)), 
#'   marker_name = paste0("marker", sprintf("%02d", 1:20)), 
#'   marker_class = factor(c(rep("type", 10), rep("state", 10)), 
#'                         levels = c("type", "state", "none")), 
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
#' design <- createDesignMatrix(experiment_info, cols_design = "group_id")
#' 
#' # Create contrast matrix
#' contrast <- createContrast(c(0, 1))
#' 
#' # Test for differential abundance (DA) of clusters
#' res_DA <- testDA_edgeR(d_counts, design, contrast)
#' 
testDA_edgeR <- function(d_counts, design, contrast, 
                         min_cells = 3, min_samples = NULL, 
                         normalize = FALSE, norm_factors = "TMM", 
                         trend_method = "none") {
  
  if (is.null(min_samples)) {
    min_samples <- ncol(d_counts) / 2
  }
  
  counts <- assays(d_counts)[["counts"]]
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
  }
  
  # note: when using DGEList object, column sums are automatically used for library sizes
  if (normalize) {
    y <- DGEList(counts, norm.factors = norm_factors)
  } else {
    y <- DGEList(counts)
  }
  
  # estimate dispersions
  # (note: using 'trend.method = "none"' by default)
  y <- estimateDisp(y, design, trend.method = trend_method)
  
  # fit models
  fit <- glmFit(y, design)
  
  # likelihood ratio tests
  lrt <- glmLRT(fit, contrast = contrast)
  
  # results
  top <- edgeR::topTags(lrt, n = Inf, adjust.method = "BH", sort.by = "none")
  
  if (!all(rownames(top) == cluster_id)) {
    stop("cluster labels do not match")
  }
  
  # return results in 'rowData' of new 'SummarizedExperiment' object
  
  # fill in any missing rows (filtered clusters) with NAs
  row_data <- as.data.frame(matrix(as.numeric(NA), nrow = nlevels(cluster_id), ncol = ncol(top)))
  colnames(row_data) <- colnames(top)
  cluster_id_nm <- as.numeric(cluster_id)
  row_data[cluster_id_nm, ] <- top$table
  
  row_data <- cbind(cluster_id = rowData(d_counts)$cluster_id, row_data)
  
  colnames(row_data)[colnames(row_data) == "PValue"] <- "p_val"
  colnames(row_data)[colnames(row_data) == "FDR"] <- "p_adj"
  
  res <- d_counts
  
  rowData(res) <- row_data
  
  # return normalization factors in 'metadata'
  if (normalize) {
    metadata(res)$norm_factors <- norm_factors
  }
  
  res
}


