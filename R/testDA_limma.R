#' Test for differential abundance: method 'diffcyt-DA-limma'
#' 
#' Calculate tests for differential abundance of clusters using method 'diffcyt-DA-limma'
#' 
#' Calculates tests for differential abundance of clusters, using empirical Bayes
#' moderation of cluster-level variances to improve power.
#' 
#' The \code{\link[limma]{limma}} package (Ritchie et al. 2015, \emph{Nucleic Acids
#' Research}) is used to fit linear models and calculate empirical Bayes moderated tests.
#' Empirical Bayes methods improve statistical power by sharing information on variability
#' (i.e. variance across samples for a single cluster) between clusters. Since count data
#' are often heteroscedastic, we use the  \code{\link[limma]{voom}} method (Law et al.
#' 2014, \emph{Genome Biology}) to transform the raw cluster cell counts and estimate
#' observation-level weights to stabilize the mean-variance relationship. Diagnostic plots
#' are shown if \code{plot = TRUE}.
#' 
#' The \code{group_IDs} argument specifies the group membership of each sample (e.g.
#' diseased vs. healthy, or treated vs. untreated), and the (optional) \code{contrast}
#' argument specifies the differential comparison of interest (e.g. group 2 vs. group 1).
#' If no \code{contrast} is provided, the default is to compare the second vs. first level
#' of \code{group_IDs}.
#' 
#' Since the \code{limma} package fits linear models for each cluster, it is possible to
#' specify flexible experimental designs. For example, the \code{group_IDs} labels may
#' contain multiple groups or conditions; in this case the \code{contrast} argument should
#' be provided to ensure that the differential tests are calculated for the contrast of
#' interest. Batch effects and continuous covariates can be provided with the
#' \code{batch_IDs} and \code{covariates} arguments; these are then added to the design
#' matrices in order to remove their effects. Paired experimental designs can be specified
#' by providing the \code{block_IDs} argument; the \code{limma}
#' \code{\link[limma]{duplicateCorrelation}} methodology is used in this case.
#' 
#' Filtering: Clusters are kept for differential testing if they have at least
#' \code{min_cells} cells in at least \code{min_samples} samples in at least one
#' condition. This removes clusters with very low cell counts across conditions, which
#' improves statistical power.
#' 
#' 
#' @param d_counts \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'   containing cluster cell counts, from \code{\link{calcCounts}}.
#' 
#' @param group_IDs Vector or factor of group membership labels for each sample (e.g.
#'   diseased vs. healthy, or treated vs. untreated). Vectors are converted to factors
#'   internally. The user must ensure that this is identical to the \code{group_IDs}
#'   provided previously to \code{\link{prepareData}}.
#' 
#' @param contrast Contrast specifying the differential comparison of interest for
#'   testing. This should be constructed using \code{\link[limma]{makeContrasts}} from the
#'   \code{\link[limma]{limma}} package. If not provided, the default is to compare the
#'   second vs. first level of \code{group_IDs}.
#' 
#' @param batch_IDs (Optional) Vector or factor of batch IDs. Batch effects are removed by
#'   adding the batch IDs to the linear model design matrices.
#' 
#' @param covariates (Optional) Numeric matrix of continuous covariates, to be added to
#'   the linear model design matrices in order to remove their effects.
#' 
#' @param block_IDs (Optional) Vector or factor of block IDs, for paired experimental
#'   designs. If provided, the \code{limma} \code{\link[limma]{duplicateCorrelation}}
#'   methodology is used to account for the paired design.
#' 
#' @param min_cells Filtering parameter. Default = 3. Clusters are kept for differential
#'   testing if they have at least \code{min_cells} cells in at least \code{min_samples}
#'   samples in at least one condition.
#' 
#' @param min_samples Filtering parameter. Default = \code{min(table(group_IDs)) - 1},
#'   i.e. one less than the number of samples in the smallest group. Clusters are kept for
#'   differential testing if they have at least \code{min_cells} cells in at least
#'   \code{min_samples} samples in at least one condition.
#' 
#' @param plot Whether to save diagnostic plots for the \code{limma}
#'   \code{\link[limma]{voom}} transformations. Default = TRUE.
#' 
#' @param path Path for diagnostic plots. Default = current working directory.
#' 
#' 
#' @return Returns a new \code{\link[SummarizedExperiment]{SummarizedExperiment}} object,
#'   with differential test results stored in the \code{rowData} slot. Results include raw
#'   p-values and adjusted p-values from the \code{limma} empirical Bayes moderated tests,
#'   which can be used to rank clusters by evidence for differential abundance. The
#'   results can be accessed with the \code{\link[SummarizedExperiment]{rowData}} accessor
#'   function.
#' 
#' 
#' @importFrom SummarizedExperiment assays rowData 'rowData<-' colData 'colData<-'
#' @importFrom limma makeContrasts contrasts.fit voom duplicateCorrelation lmFit eBayes
#'   plotSA topTable
#' @importFrom stats model.matrix
#' @importFrom methods as is
#' @importFrom grDevices pdf
#' @importFrom graphics plot
#' 
#' @export
#' 
#' @seealso to do
#' 
#' @examples
#' # to do
#' 
testDA_limma <- function(d_counts, group_IDs, contrast = NULL, 
                         batch_IDs = NULL, covariates = NULL, block_IDs = NULL, 
                         min_cells = 3, min_samples = NULL, plot = TRUE, path = ".") {
  
  if (!is.factor(group_IDs)) {
    group_IDs <- factor(group_IDs, levels = unique(group_IDs))
  }
  if (!is.null(batch_IDs) & !is.factor(batch_IDs)) {
    batch_IDs <- factor(batch_IDs, levels = unique(batch_IDs))
  }
  if (!is.null(block_IDs) & !is.factor(block_IDs)) {
    block_IDs <- factor(block_IDs, levels = unique(block_IDs))
  }
  
  if (!is.null(covariates) & !is.matrix(covariates) & !is.numeric(covariates)) {
    stop("'covariates' must be provided as a numeric matrix, with one column for each covariate")
  }
  
  if (is.null(min_samples)) {
    min_samples <- min(table(group_IDs)) - 1
  }
  
  counts <- assays(d_counts)[[1]]
  cluster <- rowData(d_counts)$cluster
  
  # filtering: keep clusters with at least 'min_cells' cells in at least 'min_samples'
  # samples in at least one condition
  ix_keep <- rep(FALSE, length(cluster))
  tf <- counts >= min_cells
  for (g in seq_along(levels(group_IDs))) {
    grp <- group_IDs == levels(group_IDs)[g]
    ix_keep[rowSums(tf[, grp, drop = FALSE]) >= min_samples] <- TRUE
  }
  
  counts <- counts[ix_keep, ]
  cluster <- cluster[ix_keep]
  
  # create design matrix
  # note: allows batch effects and covariates
  if (!is.null(batch_IDs) & !is.null(covariates)) {
    design <- model.matrix(~ 0 + group_IDs + batch_IDs + covariates)
  } else if (!is.null(batch_IDs)) {
    design <- model.matrix(~ 0 + group_IDs + batch_IDs)
  } else if (!is.null(covariates)) {
    design <- model.matrix(~ 0 + group_IDs + covariates)
  } else {
    design <- model.matrix(~ 0 + group_IDs)
  }
  
  # specify contrast of interest
  # (i.e. multiple conditions; select which one to compare to reference)
  # note: if not specified, default is to compare 2nd vs. 1st level of 'group_IDs'
  if (is.null(contrast)) {
    levs <- paste0("group_IDs", levels(group_IDs))
    my_args <- list(paste(as.character(levs[2]), "-", as.character(levs[1])), levels = design)
    contrast <- do.call(makeContrasts, my_args)
  }
  
  # voom transformation and weights
  if (plot) pdf(file.path(path, "voom_before.pdf"), width = 6, height = 6)
  v <- voom(counts, design, plot = TRUE)
  if (plot) dev.off()
  
  # estimate correlation between paired samples
  # (note: paired designs only; >2 measurements per sample not allowed)
  if (!is.null(block_IDs)) {
    dupcor <- duplicateCorrelation(v, design, block = block_IDs)
  } else {
    dupcor <- NULL
  }
  
  # fit linear models
  vfit <- lmFit(v, design, block = block_IDs, correlation = dupcor$consensus.correlation)
  vfit <- contrasts.fit(vfit, contrast)
  
  # calculate empirical Bayes moderated tests
  efit <- eBayes(vfit)
  
  if (plot) pdf(file.path(path, "voom_after.pdf"), width = 6, height = 6)
  plotSA(efit)
  if (plot) dev.off()
  
  # results
  top <- topTable(efit, coef = 1, number = Inf, adjust.method = "BH", sort.by = "none")
  
  if (!all(rownames(top) == cluster)) {
    stop("cluster labels do not match")
  }
  
  # return results in 'rowData' of new 'SummarizedExperiment' object
  
  # fill in any missing rows (filtered clusters) with NAs
  row_data <- as.data.frame(matrix(NA, nrow = nlevels(cluster), ncol = ncol(top)))
  colnames(row_data) <- colnames(top)
  cluster_nm <- as.numeric(cluster)
  row_data[cluster_nm, ] <- top
  
  row_data <- cbind(data.frame(cluster = as.numeric(levels(cluster))), row_data)
  
  # store additional sample information in 'colData'
  col_data <- cbind(colData(d_counts), data.frame(group_IDs))
  
  res <- d_counts
  
  rowData(res) <- row_data
  colData(res) <- col_data
  
  res
}


