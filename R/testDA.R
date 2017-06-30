#' Test for differential abundance
#' 
#' 
#' @param 
#' 
#' 
#' @return 
#' 
#' 
#' @importFrom SummarizedExperiment assays rowData 'rowData<-' colData 'colData<-'
#' @importFrom limma voom lmFit eBayes plotSA topTable
#' @importFrom stats model.matrix
#' @importFrom methods as is
#' @importFrom grDevices pdf
#' @importFrom graphics plot
#' 
#' @export
#' 
#' @examples
#' 
#' 
testDA_limma <- function(d_counts, group_IDs, 
                         min_cells = 5, min_samples = NULL, 
                         plot = FALSE, path = ".") {
  
  if (!is.factor(group_IDs)) group_IDs <- factor(group_IDs, levels = unique(group_IDs))
  
  if (is.null(min_samples)) {
    min_samples <- min(table(group_IDs)) - 1
  }
  
  counts <- assays(d_counts)[[1]]
  cluster <- rowData(d_counts)$cluster
  
  # filtering
  grp <- group_IDs == levels(group_IDs)[1]
  tf <- counts >= min_cells
  ix_keep <- (rowSums(tf[, grp]) >= min_samples) | (rowSums(tf[, !grp]) >= min_samples)
  
  counts <- counts[ix_keep, ]
  cluster <- cluster[ix_keep]
  
  # model matrix
  mm <- model.matrix(~ group_IDs)
  
  # voom transformation and weights
  if (plot) {
    pdf(file.path(path, "testDA_mean_var_pre_voom.pdf"), width = 6, height = 6)
    v <- voom(counts, design = mm, plot = TRUE)
    dev.off()
  } else {
    v <- voom(counts, design = mm, plot = FALSE)
  }
  
  # fit linear models
  vfit <- lmFit(v, design = mm)
  
  # calculate empirical Bayes moderated tests
  efit <- eBayes(vfit)
  
  if (plot) {
    pdf(file.path(path, "testDA_mean_var_post_voom.pdf"), width = 6, height = 6)
    plotSA(efit)
    dev.off()
  }
  
  # return new 'SummarizedExperiment' object with results stored in 'rowData'
  
  top <- topTable(efit, coef = 2, number = Inf, adjust.method = "BH", sort.by = "none")
  if (!all(rownames(top) == cluster)) {
    stop("cluster labels do not match")
  }
  
  # fill in missing rows (filtered clusters) with NAs
  row_data <- as.data.frame(matrix(NA, nrow = nlevels(cluster), ncol = ncol(top)))
  colnames(row_data) <- colnames(top)
  cluster_nm <- as.numeric(cluster)
  row_data[cluster_nm, ] <- top
  
  row_data <- cbind(data.frame(cluster = as.numeric(levels(cluster))), row_data)
  
  # also store additional sample information in 'colData'
  col_data <- cbind(colData(d_counts), data.frame(group_IDs))
  
  res_DA <- d_counts
  
  rowData(res_DA) <- row_data
  colData(res_DA) <- col_data
  
  res_DA
}


