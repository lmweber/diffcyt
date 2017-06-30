#' Test for differential abundance
#' 
#' Calculate tests for differential abundance of clusters
#' 
#' method: limma-voom
#' 
#' 
#' @param d_counts to do
#' @param group_IDs to do
#' @param contrast to do
#' @param batch_IDs to do
#' @param block_IDs to do
#' @param covariates to do
#' @param min_cells to do
#' @param min_samples to do
#' @param path to do
#' 
#' 
#' @return returns
#' 
#' 
#' @importFrom SummarizedExperiment assays rowData 'rowData<-' colData 'colData<-'
#' @importFrom limma makeContrasts voom duplicateCorrelation lmFit eBayes plotSA topTable
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
                         batch_IDs = NULL, block_IDs = NULL, covariates = NULL, 
                         min_cells = 3, min_samples = NULL, path = ".") {
  
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
  # note: allows batch effects and covariates (these can also be NULL)
  design <- model.matrix(~ 0 + group_IDs + batch_IDs + covariates)
  
  # specify contrast of interest
  # (i.e. multiple conditions; select which one to compare to reference)
  # note: if not specified, default is to compare 2nd vs. 1st level of 'group_IDs'
  if (is.null(contrast)) {
    levs <- paste0("group_IDs", levels(group_IDs))
    my_args <- list(paste(as.character(levs[2]), "-", as.character(levs[1])), levels = design)
    contrast <- do.call(makeContrasts, my_args)
  }
  
  # voom transformation and weights
  pdf(file.path(path, "voom_before.pdf"), width = 6, height = 6)
  v <- voom(counts, design, plot = TRUE)
  dev.off()
  
  # estimate correlation between paired samples
  # (note: paired designs only; >2 measurements not allowed)
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
  
  pdf("./voom_after.pdf", width = 6, height = 6)
  plotSA(efit)
  dev.off()
  
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
  
  res_DA <- d_counts
  
  rowData(res_DA) <- row_data
  colData(res_DA) <- col_data
  
  res_DA
}


