#' Test for differential functional states: method 'diffcyt-DS-limma'
#' 
#' Calculate tests for differential functional states of clusters using method
#' 'diffcyt-DS-limma'
#' 
#' Calculates tests for differential functional states of clusters (i.e. differential
#' expression of functional state markers within clusters). Clusters are defined using
#' cell type markers, and functional states are characterized by the median transformed
#' expression of functional state markers. Uses empirical Bayes moderation of
#' cluster-level variances to improve statistical power.
#' 
#' This method uses the \code{\link[limma]{limma}} package (Ritchie et al. 2015,
#' \emph{Nucleic Acids Research}) to fit linear models and calculate empirical Bayes
#' moderated tests at the cluster level. Empirical Bayes methods improve statistical power
#' by sharing information on variability (i.e. variance across samples for a single
#' cluster) between clusters. We use the option \code{trend = TRUE} in the \code{eBayes}
#' fitting function in order to stabilize the mean-variance relationship. Diagnostic plots
#' are shown if \code{plot = TRUE}.
#' 
#' The experimental design must be specified using a design matrix, which can be created
#' with \code{\link{createDesignMatrix}}. Flexible experimental designs are possible,
#' including blocking (e.g. paired designs), batch effects, and continuous covariates. See
#' \code{\link{createDesignMatrix}} for more details.
#' 
#' For paired designs, either fixed effects or random effects can be used. Fixed effects
#' are simpler, but random effects may improve power for unbalanced designs. To use fixed
#' effects, provide the block IDs (e.g. patient IDs) to \code{\link{createDesignMatrix}}.
#' To use random effects, provide the \code{block_IDs} argument here instead. This will
#' make use of the \code{limma} \code{\link[limma]{duplicateCorrelation}} methodology;
#' note that block IDs should not be included in the design matrix in this case.
#' 
#' The contrast matrix specifying the contrast of interest can be created with
#' \code{\link{createContrast}}. See \code{\link{createContrast}} for more details.
#' 
#' Filtering: Clusters are kept for differential testing if they have at least
#' \code{min_cells} cells in at least \code{min_samples} samples in at least one
#' condition. This removes clusters with very low cell counts across conditions, which
#' improves power.
#' 
#' 
#' @param d_counts \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'   containing cluster cell counts, from \code{\link{calcCounts}}.
#' 
#' @param d_medians \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'   containing cluster medians (median expression of each marker for each cluster-sample
#'   combination), from \code{\link{calcMedians}}.
#' 
#' @param design Design matrix, created with \code{\link{createDesignMatrix}}. See
#'   \code{\link{createDesignMatrix}} for details.
#' 
#' @param contrast Contrast matrix, created with \code{\link{createContrast}}. See
#'   \code{\link{createContrast}} for details.
#' 
#' @param block_IDs (Optional) Vector or factor of block IDs (e.g. patient IDs) for paired
#'   experimental designs, to be included as random effects. If provided, the block IDs
#'   will be included as random effects using the \code{limma}
#'   \code{\link[limma]{duplicateCorrelation}} methodology. Alternatively, block IDs can
#'   be included as fixed effects in the design matrix (\code{\link{createDesignMatrix}}).
#'   See details.
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
#' @param plot Whether to save diagnostic plot. Default = TRUE.
#' 
#' @param path Path for diagnostic plot. Default = current working directory.
#' 
#' 
#' @return Returns a new \code{\link[SummarizedExperiment]{SummarizedExperiment}} object,
#'   where rows = cluster-marker combinations, and columns = samples. In the rows,
#'   clusters are repeated for each functional state marker (i.e. the sheets or 'assays'
#'   from the previous \code{d_medians} object are stacked into a single matrix).
#'   Differential test results are stored in the \code{rowData} slot. Results include raw
#'   p-values and adjusted p-values from the \code{limma} empirical Bayes moderated tests,
#'   which can be used to rank cluster-marker combinations by evidence for differential
#'   functional states. The results can be accessed with the
#'   \code{\link[SummarizedExperiment]{rowData}} accessor function.
#' 
#' 
#' @importFrom SummarizedExperiment assays rowData 'rowData<-' colData 'colData<-'
#' @importFrom limma contrasts.fit duplicateCorrelation lmFit eBayes plotSA topTable
#' @importFrom methods as is
#' @importFrom grDevices pdf
#' @importFrom graphics plot
#' 
#' @export
#' 
#' @examples
#' # A full workflow example demonstrating the use of each function in the 'diffcyt'
#' # pipeline on an experimental data set is available in the package vignette.
#' 
testDS_limma <- function(d_counts, d_medians, design, contrast, 
                         block_IDs = NULL, 
                         min_cells = 3, min_samples = NULL, 
                         plot = TRUE, path = ".") {
  
  if (!is.null(block_IDs) & !is.factor(block_IDs)) {
    block_IDs <- factor(block_IDs, levels = unique(block_IDs))
  }
  
  group_IDs <- colData(d_counts)$group_IDs
  
  if (is.null(min_samples)) {
    min_samples <- min(table(group_IDs)) - 1
  }
  
  # vector identifying functional state markers
  id_state_markers <- metadata(d_medians)$id_state_markers
  
  # note: counts are only required for filtering
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
  
  cluster <- cluster[ix_keep]
  
  # extract medians and create concatenated matrix
  state_names <- names(assays(d_medians))[id_state_markers]
  meds <- do.call("rbind", {
    lapply(as.list(assays(d_medians)[state_names]), function(a) a[cluster, ])
  })
  
  meds_all <- do.call("rbind", as.list(assays(d_medians)[state_names]))
  
  # limma pipeline
  
  # estimate correlation between paired samples
  # (note: paired designs only; >2 measures per sample not allowed)
  if (!is.null(block_IDs)) {
    dupcor <- duplicateCorrelation(meds, design, block = block_IDs)
  }
  
  # fit linear models
  if (!is.null(block_IDs)) {
    message("Fitting linear models with random effects term for 'block_IDs'.")
    fit <- lmFit(meds, design, block = block_IDs, correlation = dupcor$consensus.correlation)
  } else {
    fit <- lmFit(meds, design)
  }
  fit <- contrasts.fit(fit, contrast)
  
  # calculate empirical Bayes moderated tests (note: using 'trend = TRUE' for mean-variance relationship)
  efit <- eBayes(fit, trend = TRUE)
  
  if (plot) pdf(file.path(path, "SA_plot.pdf"), width = 6, height = 6)
  plotSA(efit)
  if (plot) dev.off()
  
  # results
  top <- topTable(efit, coef = 1, number = Inf, adjust.method = "BH", sort.by = "none")
  
  if (!all(top$ID %in% cluster)) {
    stop("cluster labels do not match")
  }
  
  # return results in 'rowData' of new 'SummarizedExperiment' object
  
  # fill in any missing rows (filtered clusters) with NAs
  row_data <- as.data.frame(matrix(as.numeric(NA), 
                                   nrow = nlevels(cluster) * length(state_names), 
                                   ncol = ncol(top)))
  colnames(row_data) <- colnames(top)
  
  cluster_nm <- as.numeric(cluster)
  s <- seq(0, nlevels(cluster) * (length(state_names) - 1), by = nlevels(cluster))
  r1 <- rep(cluster_nm, length(state_names))
  r2 <- rep(s, each = length(cluster_nm))
  
  stopifnot(length(s) == length(state_names))
  stopifnot(length(r1) == length(r2))
  
  rows <- r1 + r2
  row_data[rows, ] <- top
  
  # include cluster IDs and marker names
  clus <- factor(rep(levels(cluster), length(state_names)), levels = levels(cluster))
  stat <- factor(rep(state_names, each = length(levels(cluster))), levels = state_names)
  stopifnot(length(clus) == nrow(row_data), length(stat) == nrow(row_data))
  
  row_data <- cbind(data.frame(cluster = clus, marker = stat), row_data)
  
  col_data <- colData(d_medians)
  
  # return object
  res <- SummarizedExperiment(
    meds_all, 
    rowData = row_data, 
    colData = col_data
  )
  
  res
}


