#' Plot heatmap
#' 
#' Plot heatmap showing top clusters or cluster-marker combinations
#' 
#' Display heatmap to visualize results for the top (most highly significant) detected
#' clusters or cluster-marker combinations.
#' 
#' For DA tests, the heatmap consists of the following panels:
#' 
#' \itemize{
#' \item median (arcsinh-transformed) expression (across all samples) for 'cell type'
#' markers
#' \item cluster abundances by sample
#' \item row annotation indicating significant detected clusters
#' }
#' 
#' For DS tests, the heatmap consists of:
#' 
#' \itemize{
#' \item median (arcsinh-transformed) expression (across all samples) for 'cell type'
#' markers
#' \item median (arcsinh-transformed) expression (across all samples) for 'cell state'
#' markers
#' \item median (arcsinh-transformed) expression (by sample) for 'cell state' markers for
#' the top cluster-marker combinations
#' \item row annotation indicating significant detected cluster-marker combinations
#' }
#' 
#' Heatmaps are generated using the \code{ComplexHeatmap} package (Gu et al., 2016), and
#' color scales are generated using the \code{circlize} package (Gu et al., 2014). Both
#' packages are available from Bioconductor.
#' 
#' 
#' @param out Output object from \code{\link{diffcyt}} wrapper function, containing
#'   results object \code{res} and data objects \code{d_se}, \code{d_counts},
#'   \code{d_medians}, and \code{d_medians_by_cluster_marker}. Alternatively, the results
#'   and data objects can be provided individually.
#' 
#' @param analysis_type Whether to plot heatmap for differential abundance (DA) or
#'   differential state (DS) test results.
#' 
#' @param top_n Number of top clusters (DA tests) or cluster-marker combinations (DS
#'   tests) to display. Default = 20.
#' 
#' @param threshold Threshold for significant adjusted p-values. Default = 0.1.
#' 
#' @param res Object containing differential test results. Alternatively, the combined
#'   output object from the wrapper function \code{\link{diffcyt}} can be provided.
#' 
#' @param d_se Data object. Alternatively, the combined output object from the wrapper
#'   function \code{\link{diffcyt}} can be provided.
#' 
#' @param d_counts Data object. Alternatively, the combined output object from the wrapper
#'   function \code{\link{diffcyt}} can be provided.
#' 
#' @param d_medians Data object. (Required for DS tests only.) Alternatively, the combined
#'   output object from the wrapper function \code{\link{diffcyt}} can be provided.
#' 
#' @param d_medians_by_cluster_marker Data object. Alternatively, the combined output
#'   object from the wrapper function \code{\link{diffcyt}} can be provided.
#' 
#' @param sample_order (Optional) Custom ordering for samples (columns) in right-hand
#'   panel of heatmap. (This is useful when the default ordering does not group samples by
#'   condition; e.g. samples are ordered alphabetically by sample IDs instead.)
#' 
#' 
#' @return Displays a heatmap.
#' 
#' 
#' @importFrom ComplexHeatmap Heatmap columnAnnotation rowAnnotation anno_text
#'   max_text_width draw '+.AdditiveUnit'
#' @importFrom circlize colorRamp2
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom S4Vectors metadata
#' @importFrom stats quantile
#' @importFrom grid unit unit.c gpar
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
#' # Add differential states (DS) signal
#' ix_DS <- 901:1000
#' ix_cols_DS <- 19:20
#' d_input[[1]][ix_DS, ix_cols_type] <- d_random(n = 1000, mean = 3, ncol = 10)
#' d_input[[2]][ix_DS, ix_cols_type] <- d_random(n = 1000, mean = 3, ncol = 10)
#' d_input[[3]][ix_DS, c(ix_cols_type, ix_cols_DS)] <- d_random(n = 1200, mean = 3, ncol = 12)
#' d_input[[4]][ix_DS, c(ix_cols_type, ix_cols_DS)] <- d_random(n = 1200, mean = 3, ncol = 12)
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
#' # Create design matrix
#' design <- createDesignMatrix(experiment_info, cols_design = 2)
#' 
#' # Create contrast matrix
#' contrast <- createContrast(c(0, 1))
#' 
#' # Test for differential abundance (DA) of clusters (using default method 'diffcyt-DA-edgeR')
#' out_DA <- diffcyt(d_input, experiment_info, marker_info, 
#'                   design = design, contrast = contrast, 
#'                   analysis_type = "DA", method_DA = "diffcyt-DA-edgeR", 
#'                   seed_clustering = 123, verbose = FALSE)
#' 
#' # Test for differential states (DS) within clusters (using default method 'diffcyt-DS-limma')
#' out_DS <- diffcyt(d_input, experiment_info, marker_info, 
#'                   design = design, contrast = contrast, 
#'                   analysis_type = "DS", method_DS = "diffcyt-DS-limma", 
#'                   seed_clustering = 123, plot = FALSE, verbose = FALSE)
#' 
#' # Display results for top DA clusters
#' topTable(out_DA, format_vals = TRUE)
#' 
#' # Display results for top DS cluster-marker combinations
#' topTable(out_DS, format_vals = TRUE)
#' 
#' # Plot heatmap for DA tests
#' plotHeatmap(out_DA, analysis_type = "DA")
#' 
#' # Plot heatmap for DS tests
#' plotHeatmap(out_DS, analysis_type = "DS")
#' 
plotHeatmap <- function(out = NULL, analysis_type = c("DA", "DS"), top_n = 20, threshold = 0.1, 
                        res = NULL, d_se = NULL, d_counts = NULL, d_medians = NULL, d_medians_by_cluster_marker = NULL, 
                        sample_order = NULL) {
  
  if (!is.null(out)) {
    res <- out$res
    d_se <- out$d_se
    d_counts <- out$d_counts
    d_medians <- out$d_medians
    d_medians_by_cluster_marker <- out$d_medians_by_cluster_marker
  }
  
  analysis_type <- match.arg(analysis_type, choices = c("DA", "DS"))
  
  
  is_marker <- colData(d_medians_by_cluster_marker)$marker_class != "none"
  is_celltype_marker <- colData(d_medians_by_cluster_marker)$marker_class == "type"
  is_state_marker <- colData(d_medians_by_cluster_marker)$marker_class == "state"
  
  
  # ------------------------------------------------------------------------------
  # heatmap: main panel (DA tests and DS tests): expression of 'cell type' markers
  # ------------------------------------------------------------------------------
  
  d_heatmap <- 
    assay(d_medians_by_cluster_marker)[, is_marker, drop = FALSE]
  
  d_heatmap_celltype <- 
    assay(d_medians_by_cluster_marker)[, is_celltype_marker, drop = FALSE]
  
  # arrange alphabetically
  d_heatmap_celltype <- d_heatmap_celltype[, order(colnames(d_heatmap_celltype)), drop = FALSE]
  
  if (analysis_type == "DA") {
    stopifnot(nrow(d_heatmap) == nrow(rowData(res)$cluster_id), 
              nrow(d_heatmap_celltype) == nrow(rowData(res)$cluster_id), 
              all(rownames(d_heatmap) == rowData(res)$cluster_id), 
              all(rownames(d_heatmap_celltype) == rowData(res)$cluster_id))
  } else if (analysis_type == "DS") {
    stopifnot(nrow(d_heatmap) == nlevels(rowData(res)$cluster_id), 
              nrow(d_heatmap_celltype) == nlevels(rowData(res)$cluster_id), 
              all(rownames(d_heatmap) %in% rowData(res)$cluster_id), 
              all(rownames(d_heatmap_celltype) %in% rowData(res)$cluster_id))
  }
  
  # results for top clusters
  d_top <- topTable(res, top_n = top_n, format_vals = FALSE)
  
  # top 'n' clusters
  d_heatmap_celltype <- d_heatmap_celltype[match(d_top$cluster_id, rownames(d_heatmap_celltype)), , drop = FALSE]
  
  stopifnot(nrow(d_heatmap_celltype) == nrow(d_top), 
            all(rownames(d_heatmap_celltype) == d_top$cluster_id))
  
  # color scale: 1%, 50%, 99% percentiles across all medians and all markers
  colors <- colorRamp2(
    quantile(assay(d_medians_by_cluster_marker)[, is_marker], 
             c(0.01, 0.5, 0.99), na.rm = TRUE), 
    c("royalblue3", "white", "tomato2")
  )
  
  # note: no additional scaling (using asinh-transformed values directly)
  ht_main <- Heatmap(
    d_heatmap_celltype, col = colors, name = "expression", 
    row_title = "clusters", row_title_gp = gpar(fontsize = 14), 
    column_title = "markers (cell type)", column_title_side = "bottom", column_title_gp = gpar(fontsize = 14), 
    column_names_gp = gpar(fontsize = 12), 
    heatmap_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12)), 
    cluster_columns = FALSE, row_names_side = "left", row_names_gp = gpar(fontsize = 11), 
    clustering_distance_rows = "euclidean", clustering_method_rows = "median"
  )
  
  
  # --------------------------------------------------------------
  # heatmap: second panel (DA tests): cluster abundances by sample
  # --------------------------------------------------------------
  
  if (analysis_type == "DA") {
    
    stopifnot(nrow(d_counts) == nrow(rowData(res)), 
              all(rownames(assay(d_counts)) == rownames(rowData(res))))
    
    d_abundance <- assay(d_counts)[d_top$cluster_id, , drop = FALSE]
    
    stopifnot(nrow(d_abundance) == nrow(d_heatmap_celltype), 
              all(rownames(d_abundance) == rownames(d_heatmap_celltype)))
    
    # color scale: full range
    colors_counts <- colorRamp2(range(d_abundance), c("#132a13", "yellow"))
    
    # note: row ordering is automatically matched when multiple heatmaps are combined
    ht_abundance <- Heatmap(
      d_abundance, col = colors_counts, name = "n_cells", 
      column_title = "samples", column_title_side = "bottom", column_title_gp = gpar(fontsize = 14), 
      column_names_gp = gpar(fontsize = 12), 
      heatmap_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12)), 
      column_order = sample_order, cluster_columns = FALSE, 
      show_row_names = FALSE
    )
  }
  
  
  # ------------------------------------------------------------------------------
  # heatmap: second panel (DS tests): expression of 'cell state' markers by sample
  # ------------------------------------------------------------------------------
  
  if (analysis_type == "DS") {
    
    # create data frame of expression values of top 'n' cluster-marker combinations by sample
    top_clusters <- as.list(d_top$cluster_id)
    top_markers <- as.character(d_top$marker_id)
    
    assays_ordered <- assays(d_medians)[top_markers]
    
    d_markers <- mapply(function(a, cl) {
      a[cl, , drop = FALSE]
    }, assays_ordered, top_clusters, SIMPLIFY = FALSE)
    
    d_markers <- do.call("rbind", d_markers)
    
    stopifnot(all(rownames(d_markers) == rownames(d_top)), 
              all(rownames(d_markers) == rownames(d_heatmap_celltype)), 
              nrow(d_markers) == nrow(d_heatmap_celltype), 
              all(top_markers == d_top$marker))
    
    rownames(d_markers) <- paste0(top_markers, " (", rownames(d_markers), ")")
    
    # color scale: full range
    colors_markers <- colorRamp2(range(d_markers, na.rm = TRUE), c("navy", "yellow"))
    
    # note: row ordering is automatically matched when multiple heatmaps are combined
    ht_markers <- Heatmap(
      d_markers, col = colors_markers, name = "expression\n(by sample)",
      column_title = "samples", column_title_side = "bottom", column_title_gp = gpar(fontsize = 14), 
      column_names_gp = gpar(fontsize = 12), 
      heatmap_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12)), 
      column_order = sample_order, cluster_columns = FALSE, 
      show_row_names = TRUE, row_names_side = "right", row_names_gp = gpar(fontsize = 11)
    )
  }
  
  
  # ---------------------------------
  # row annotation: adjusted p-values
  # ---------------------------------
  
  # identify column of adjusted p-values
  ix_p_adj <- which(colnames(d_top) == "p_adj")
  
  # significant differential clusters or cluster-marker combinations
  sig <- d_top[, ix_p_adj] <= threshold
  # set filtered clusters or cluster-marker combinations to FALSE
  sig[is.na(sig)] <- FALSE
  
  # set up data frame
  if (analysis_type == "DA") {
    d_sig <- data.frame(cluster_id = d_top$cluster_id, 
                        sig = as.numeric(sig))
  } else if (analysis_type == "DS") {
    d_sig <- data.frame(cluster_id = d_top$cluster_id, 
                        marker = d_top$marker_id, 
                        sig = as.numeric(sig))
  }
  
  stopifnot(nrow(d_sig) == nrow(d_top), 
            nrow(d_sig) == nrow(d_heatmap_celltype))
  
  # add row annotation
  row_annot <- data.frame(
    "significant" = factor(d_sig$sig, levels = c(0, 1), labels = c("no", "yes")), 
    check.names = FALSE
  )
  
  ha_row <- rowAnnotation(
    df = row_annot, 
    col = list("significant" = c("no" = "gray90", "yes" = "red")), 
    annotation_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12)), 
    width = unit(0.75, "cm")
  )
  
  
  # ----------------------
  # combine heatmap panels
  # ----------------------
  
  # title
  if (analysis_type == "DA") {
    ht_title <- "Results: top DA clusters"
  } else if (analysis_type == "DS") {
    ht_title <- "Results: top DS cluster-marker combinations"
  }
  
  # combine elements of heatmap
  if (analysis_type == "DA") {
    draw("+.AdditiveUnit"("+.AdditiveUnit"(ht_main, ht_abundance), ha_row), 
         column_title = ht_title, column_title_gp = gpar(fontface = "bold", fontsize = 12))
    
  } else if (analysis_type == "DS") {
    draw("+.AdditiveUnit"("+.AdditiveUnit"(ht_main, ht_markers), ha_row), 
         column_title = ht_title, column_title_gp = gpar(fontface = "bold", fontsize = 12))
  }
  
}


