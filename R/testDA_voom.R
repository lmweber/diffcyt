#' Test for differential abundance: method 'diffcyt-DA-voom'
#' 
#' Calculate tests for differential abundance of cell populations using method
#' 'diffcyt-DA-voom'
#' 
#' Calculates tests for differential abundance of clusters, using functions from the
#' \code{\link{limma}} package and \code{\link{voom}} method.
#' 
#' This method uses the \code{\link{limma}} package (Ritchie et al. 2015, \emph{Nucleic
#' Acids Research}) to fit models and calculate moderated tests at the cluster level.
#' Moderated tests improve statistical power by sharing information on variability (i.e.
#' variance across samples for a single cluster) between clusters. Since count data are
#' often heteroscedastic, we use the \code{\link{voom}} method (Law et al. 2014,
#' \emph{Genome Biology}) to transform the raw cluster cell counts and estimate
#' observation-level weights to stabilize the mean-variance relationship. Diagnostic plots
#' are shown if \code{plot = TRUE}.
#' 
#' The experimental design must be specified using a design matrix, which can be created
#' with \code{\link{createDesignMatrix}}. Flexible experimental designs are possible,
#' including blocking (e.g. paired designs), batch effects, and continuous covariates. See
#' \code{\link{createDesignMatrix}} for more details.
#' 
#' For paired designs, either fixed effects or random effects can be used. Fixed effects
#' are simpler, but random effects may improve power in data sets with unbalanced designs
#' or very large numbers of samples. To use fixed effects, provide the block IDs (e.g.
#' patient IDs) to \code{\link{createDesignMatrix}}. To use random effects, provide the
#' \code{block_id} argument here instead. This will make use of the \code{limma}
#' \code{\link{duplicateCorrelation}} methodology. Note that >2 measures per sample are
#' not possible in this case (fixed effects should be used instead). Block IDs should not
#' be included in the design matrix if the \code{limma} \code{duplicateCorrelation}
#' methodology is used.
#' 
#' The contrast matrix specifying the contrast of interest can be created with
#' \code{\link{createContrast}}. See \code{\link{createContrast}} for more details.
#' 
#' Filtering: Clusters are kept for differential testing if they have at least
#' \code{min_cells} cells in at least \code{min_samples} samples. This removes clusters
#' with very low cell counts across conditions, to improve power.
#' 
#' Normalization for the total number of cells per sample (library sizes) and total number
#' of cells per cluster is automatically performed by the \code{limma} and \code{voom}
#' functions. Optional normalization factors can also be included to adjust for
#' composition effects in the total number of counts per sample (library sizes). For
#' example, if one cell population is more abundant in one condition, while all other
#' populations have equal numbers of cells across conditions, then this effectively
#' reduces the \emph{relative} abundance of the unchanged populations in the first
#' condition, creating false positive differential abundance signals for these
#' populations. Normalization factors can be provided directly as a vector of values
#' representing relative total abundances per sample (where values >1 indicate extra cells
#' in a sample; note this is the inverse of normalization factors as defined in the
#' \code{edgeR} package). Alternatively, normalization factors can be calculated
#' automatically using the 'trimmed mean of M-values' (TMM) method from the \code{edgeR}
#' package (Robinson and Oshlack, 2010). The TMM method assumes that most populations are
#' not differentially abundant; see the \code{edgeR} User's Guide for more details.
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
#' @param block_id (Optional) Vector or factor of block IDs (e.g. patient IDs) for paired
#'   experimental designs, to be included as random effects. If provided, the block IDs
#'   will be included as random effects using the \code{limma}
#'   \code{\link{duplicateCorrelation}} methodology. Alternatively, block IDs can be
#'   included as fixed effects in the design matrix (\code{\link{createDesignMatrix}}).
#'   See details.
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
#' @param plot Whether to save diagnostic plots for the \code{limma} \code{\link{voom}}
#'   transformations. Default = TRUE.
#' 
#' @param path Path for diagnostic plots. Default = current working directory.
#' 
#' 
#' @return Returns a new \code{\link{SummarizedExperiment}} object, with differential test
#'   results stored in the \code{rowData} slot. Results include raw p-values
#'   (\code{p_val}) and adjusted p-values (\code{p_adj}) from the \code{limma} moderated
#'   tests, which can be used to rank clusters by evidence for differential abundance.
#'   Additional output columns from the \code{limma} tests are also included. The results
#'   can be accessed with the \code{\link{rowData}} accessor function.
#' 
#' 
#' @importFrom SummarizedExperiment assay rowData 'rowData<-' colData 'colData<-'
#' @importFrom limma contrasts.fit voom duplicateCorrelation lmFit eBayes plotSA topTable
#' @importFrom edgeR calcNormFactors DGEList
#' @importFrom methods as is
#' @importFrom grDevices pdf
#' @importFrom graphics plot
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
#' design <- createDesignMatrix(experiment_info, cols_design = 2)
#' 
#' # Create contrast matrix
#' contrast <- createContrast(c(0, 1))
#' 
#' # Test for differential abundance (DA) of clusters
#' res_DA <- testDA_voom(d_counts, design, contrast, plot = FALSE)
#' 
testDA_voom <- function(d_counts, design, contrast, block_id = NULL, 
                        min_cells = 3, min_samples = NULL, 
                        normalize = FALSE, norm_factors = "TMM", 
                        plot = TRUE, path = ".") {
  
  if (!is.null(block_id) & !is.factor(block_id)) {
    block_id <- factor(block_id, levels = unique(block_id))
  }
  
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
  
  # limma-voom pipeline
  
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
  
  # voom transformation and weights
  if (plot) pdf(file.path(path, "voom_before.pdf"), width = 6, height = 6)
  v <- voom(y, design, plot = TRUE)
  if (plot) dev.off()
  
  # estimate correlation between paired samples
  # (note: paired designs only; >2 measures per sample not allowed)
  if (!is.null(block_id)) {
    dupcor <- duplicateCorrelation(v, design, block = block_id)
  }
  
  # fit models
  if (!is.null(block_id)) {
    message("Fitting linear models with random effects term for 'block_id'.")
    vfit <- lmFit(v, design, block = block_id, correlation = dupcor$consensus.correlation)
  } else {
    vfit <- lmFit(v, design)
  }
  vfit <- contrasts.fit(vfit, contrast)
  
  # calculate moderated tests
  efit <- eBayes(vfit)
  
  if (plot) pdf(file.path(path, "voom_after.pdf"), width = 6, height = 6)
  plotSA(efit)
  if (plot) dev.off()
  
  # results
  top <- topTable(efit, coef = 1, number = Inf, adjust.method = "BH", sort.by = "none")
  
  if (!all(rownames(top) == cluster_id)) {
    stop("cluster labels do not match")
  }
  
  # return results in 'rowData' of new 'SummarizedExperiment' object
  
  # fill in any missing rows (filtered clusters) with NAs
  row_data <- as.data.frame(matrix(as.numeric(NA), nrow = nlevels(cluster_id), ncol = ncol(top)))
  colnames(row_data) <- colnames(top)
  cluster_id_nm <- as.numeric(cluster_id)
  row_data[cluster_id_nm, ] <- top
  
  row_data <- cbind(cluster_id = rowData(d_counts)$cluster_id, row_data)
  
  colnames(row_data)[colnames(row_data) == "P.Value"] <- "p_val"
  colnames(row_data)[colnames(row_data) == "adj.P.Val"] <- "p_adj"
  
  res <- d_counts
  
  rowData(res) <- row_data
  
  res
}


