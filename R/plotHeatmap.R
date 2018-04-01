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
#'   \code{d_medians}, and \code{d_medians_all_samples}. Alternatively, the results and
#'   data objects can be provided individually.
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
#' @param d_medians_all_samples Data object. Alternatively, the combined output object
#'   from the wrapper function \code{\link{diffcyt}} can be provided.
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
#'   sample = factor(paste0("sample", 1:4)), 
#'   group = factor(c("group1", "group1", "group2", "group2")), 
#'   stringsAsFactors = FALSE
#' )
#' 
#' marker_info <- data.frame(
#'   marker_name = paste0("marker", 1:20), 
#'   is_marker = rep(TRUE, 20), 
#'   is_type_marker = c(rep(TRUE, 10), rep(FALSE, 10)), 
#'   is_state_marker = c(rep(FALSE, 10), rep(TRUE, 10)), 
#'   stringsAsFactors = FALSE
#' )
#' 
#' # Create design matrix
#' design <- createDesignMatrix(sample_info, cols_include = 2)
#' # Create contrast matrix
#' contrast <- createContrast(c(0, 1))
#' 
#' # Test for differential abundance (DA) of clusters (using default method 'diffcyt-DA-edgeR')
#' out_DA <- diffcyt(d_input, sample_info, marker_info, 
#'                   design = design, contrast = contrast, 
#'                   analysis_type = "DA", method_DA = "diffcyt-DA-edgeR")
#' 
#' # Test for differential states (DS) within clusters (using default method 'diffcyt-DS-limma')
#' out_DS <- diffcyt(d_input, sample_info, marker_info, 
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
plotHeatmap <- function(out = NULL, analysis_type = c("DA", "DS"), top_n = 20, threshold = 0.1, 
                        res = NULL, d_se = NULL, d_counts = NULL, d_medians = NULL, d_medians_all_samples = NULL) {
  
  if (!is.null(out)) {
    res <- out$res
    d_se <- out$d_se
    d_counts <- out$d_counts
    d_medians <- out$d_medians
    d_medians_all_samples <- out$d_medians_all_samples
  }
  
  analysis_type <- match.arg(analysis_type, choices = c("DA", "DS"))
  
  
  # ------------------------------------------------------------------------------
  # heatmap: main panel (DA tests and DS tests): expression of 'cell type' markers
  # ------------------------------------------------------------------------------
  
  d_heatmap <- assay(d_medians_all_samples)[, colData(d_medians_all_samples)$is_marker, drop = FALSE]
  
  d_heatmap_celltype <- assay(d_medians_all_samples)[, colData(d_medians_all_samples)$is_type_marker, drop = FALSE]
  # arrange alphabetically
  d_heatmap_celltype <- d_heatmap_celltype[, order(colnames(d_heatmap_celltype)), drop = FALSE]
  
  if (analysis_type == "DA") {
    stopifnot(nrow(d_heatmap) == nrow(rowData(res)$cluster), 
              nrow(d_heatmap_celltype) == nrow(rowData(res)$cluster), 
              all(rownames(d_heatmap) == rowData(res)$cluster), 
              all(rownames(d_heatmap_celltype) == rowData(res)$cluster))
  } else if (analysis_type == "DS") {
    stopifnot(nrow(d_heatmap) == nlevels(rowData(res)$cluster), 
              nrow(d_heatmap_celltype) == nlevels(rowData(res)$cluster), 
              all(rownames(d_heatmap) %in% rowData(res)$cluster), 
              all(rownames(d_heatmap_celltype) %in% rowData(res)$cluster))
  }
  
  # results for top clusters
  d_top <- topClusters(res, order = TRUE, all = FALSE, top_n = top_n)
  
  # top 'n' clusters
  d_heatmap_celltype <- d_heatmap_celltype[match(d_top$cluster, rownames(d_heatmap_celltype)), , drop = FALSE]
  
  stopifnot(nrow(d_heatmap_celltype) == nrow(d_top), 
            all(rownames(d_heatmap_celltype) == d_top$cluster))
  
  # color scale: 1%, 50%, 99% percentiles across all medians and all markers
  colors <- colorRamp2(
    quantile(assay(d_medians_all_samples)[, colData(d_medians_all_samples)$is_marker], c(0.01, 0.5, 0.99), na.rm = TRUE), 
    c("royalblue3", "white", "tomato2")
  )
  
  # note: no additional scaling (using asinh-transformed values directly)
  if (analysis_type == "DA") {
    ht_main <- Heatmap(
      d_heatmap_celltype, col = colors, name = "expression", 
      row_title = "clusters", row_title_gp = gpar(fontsize = 14), 
      column_title = "markers", column_title_side = "bottom", column_title_gp = gpar(fontsize = 14), 
      column_names_gp = gpar(fontsize = 12), 
      heatmap_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12)), 
      cluster_columns = FALSE, row_names_side = "left", row_names_gp = gpar(fontsize = 11), 
      clustering_distance_rows = "euclidean", clustering_method_rows = "median"
    )
  } else if (analysis_type == "DS") {
    
    # column annotation for cell type markers
    n_celltype <- sum(metadata(d_medians_all_samples)$id_type_markers)
    
    col_annot_celltype <- data.frame(
      "marker type" = factor(c(rep("cell type", n_celltype), rep("cell state", 0)), levels = c("cell type", "cell state")), 
      check.names = FALSE
    )
    
    ha_col_celltype <- columnAnnotation(
      df = col_annot_celltype, show_legend = FALSE, 
      col = list("marker type" = c("cell type" = "gold", "cell state" = "darkgreen")), 
      colname = anno_text(colnames(d_heatmap_celltype), rot = 90, just = "right", offset = unit(1, "npc") - unit(2, "mm")), 
      annotation_height = unit.c(unit(5, "mm"), max_text_width(colnames(d_heatmap_celltype)) + unit(2, "mm")), 
      annotation_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12))
    )
    
    ht_main <- Heatmap(
      d_heatmap_celltype, col = colors, name = "expression", 
      row_title = "clusters", row_title_gp = gpar(fontsize = 14), 
      column_title = "markers (cell type)", column_title_side = "bottom", column_title_gp = gpar(fontsize = 14), 
      show_column_names = FALSE, 
      heatmap_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12)), 
      cluster_columns = FALSE, row_names_side = "left", row_names_gp = gpar(fontsize = 11), 
      clustering_distance_rows = "euclidean", clustering_method_rows = "median", 
      bottom_annotation = ha_col_celltype
    )
  }
  
  
  # --------------------------------------------------------------
  # heatmap: second panel (DA tests): cluster abundances by sample
  # --------------------------------------------------------------
  
  if (analysis_type == "DA") {
    
    stopifnot(nrow(d_counts) == nrow(rowData(res)), 
              all(rownames(assay(d_counts)) == rownames(rowData(res))))
    
    d_abundance <- assay(d_counts)[d_top$cluster, , drop = FALSE]
    
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
      cluster_columns = FALSE, show_row_names = FALSE
    )
  }
  
  
  # --------------------------------------------------------------------
  # heatmap: second panel (DS tests): expression of 'cell state' markers
  # --------------------------------------------------------------------
  
  if (analysis_type == "DS") {
    
    d_heatmap_state <- assay(d_medians_all_samples)[, colData(d_medians_all_samples)$is_state_marker, drop = FALSE]
    # arrange alphabetically
    d_heatmap_state <- d_heatmap_state[, order(colnames(d_heatmap_state)), drop = FALSE]
    
    stopifnot(nrow(d_heatmap_state) == nlevels(rowData(res)$cluster), 
              all(rownames(d_heatmap_state) %in% rowData(res)$cluster))
    
    # results for top cluster-marker combinations
    d_top <- topClusters(res, order = TRUE, all = FALSE, top_n = top_n)
    
    # top 'n' clusters
    d_heatmap_state <- d_heatmap_state[match(d_top$cluster, rownames(d_heatmap_state)), , drop = FALSE]
    
    stopifnot(nrow(d_heatmap_state) == nrow(d_top), 
              all(rownames(d_heatmap_state) == d_top$cluster))
    
    # column annotation for cell state markers
    n_state <- sum(metadata(d_medians_all_samples)$id_state_markers)
    
    col_annot_state <- data.frame(
      "marker type" = factor(c(rep("cell type", 0), rep("cell state", n_state)), levels = c("cell type", "cell state")), 
      check.names = FALSE
    )
    
    ha_col_state <- columnAnnotation(
      df = col_annot_state, show_legend = FALSE, 
      col = list("marker type" = c("cell type" = "gold", "cell state" = "darkgreen")), 
      colname = anno_text(colnames(d_heatmap_state), rot = 90, just = "right", offset = unit(1, "npc") - unit(2, "mm")), 
      annotation_height = unit.c(unit(5, "mm"), max_text_width(colnames(d_heatmap_state)) + unit(2, "mm")), 
      annotation_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12))
    )
    
    ht_state <- Heatmap(
      d_heatmap_state, col = colors, 
      show_heatmap_legend = FALSE, 
      column_title = "markers (cell state)", column_title_side = "bottom", column_title_gp = gpar(fontsize = 14), 
      show_column_names = FALSE, 
      column_names_gp = gpar(fontsize = 12), 
      cluster_columns = FALSE, show_row_names = FALSE, 
      bottom_annotation = ha_col_state
    )
  }
  
  
  # -----------------------------------------------------------------------------
  # heatmap: third panel (DS tests): expression of 'cell state' markers by sample
  # -----------------------------------------------------------------------------
  
  if (analysis_type == "DS") {
    
    # create data frame of expression values of top 'n' cluster-marker combinations by sample
    top_clusters <- as.list(d_top$cluster)
    top_markers <- as.character(d_top$marker)
    
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
    
    # column annotation (empty: need for formatting purposes)
    col_annot_markers <- data.frame(
      "state marker" = factor(rep("markers", ncol(d_markers)), levels = "markers"), 
      check.names = FALSE
    )
    
    ha_col_markers <- columnAnnotation(
      df = col_annot_markers, show_legend = FALSE, 
      col = list("state marker" = c("markers" = "gray")), 
      colname = anno_text(colnames(d_markers), rot = 90, just = "right", offset = unit(1, "npc") - unit(2, "mm")), 
      # use zero height
      annotation_height = unit.c(unit(0, "mm"), max_text_width(colnames(d_markers)) + unit(2, "mm")), 
      annotation_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12))
    )
    
    # color scale: full range
    colors_markers <- colorRamp2(range(d_markers, na.rm = TRUE), c("navy", "yellow"))
    
    # note: row ordering is automatically matched when multiple heatmaps are combined
    ht_markers <- Heatmap(
      d_markers, col = colors_markers, name = "expression\n(by sample)",
      column_title = "samples", column_title_side = "bottom", column_title_gp = gpar(fontsize = 14), 
      show_column_names = FALSE, 
      column_names_gp = gpar(fontsize = 12), 
      heatmap_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12)), 
      cluster_columns = FALSE, 
      show_row_names = TRUE, row_names_side = "right", row_names_gp = gpar(fontsize = 11), 
      bottom_annotation = ha_col_markers
    )
  }
  
  
  # ---------------------------------
  # row annotation: adjusted p-values
  # ---------------------------------
  
  # identify column of adjusted p-values
  ix_p_adj <- which(colnames(d_top) %in% c("FDR", "adj.P.Val", "p_adj"))
  
  # significant differential clusters or cluster-marker combinations
  sig <- d_top[, ix_p_adj] <= threshold
  # set filtered clusters or cluster-marker combinations to FALSE
  sig[is.na(sig)] <- FALSE
  
  # set up data frame
  if (analysis_type == "DA") {
    d_sig <- data.frame(cluster = d_top$cluster, 
                        sig = as.numeric(sig))
  } else if (analysis_type == "DS") {
    d_sig <- data.frame(cluster = d_top$cluster, 
                        marker = d_top$marker, 
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
    draw("+.AdditiveUnit"("+.AdditiveUnit"("+.AdditiveUnit"(ht_main, ht_state), ht_markers), ha_row), 
         column_title = ht_title, column_title_gp = gpar(fontface = "bold", fontsize = 12))
  }
  
}


