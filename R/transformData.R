#' Transform data
#' 
#' Transform data prior to clustering
#' 
#' Flow and mass cytometry data should be transformed prior to clustering. This function 
#' implements an 'arcsinh' transform with adjustable 'cofactor' parameter. Recommended 
#' values for the cofactor are 5 (mass cytometry, CyTOF) or 150 (fluorescence flow 
#' cytometry); see Bendall et al. (2011), \emph{Science}, Supplementary Figure S2.
#' 
#' The transform will be applied to protein marker columns only. The 
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} object created in the previous
#' step (\code{\link{prepareData}}) is assumed to contain a vector of logical entries
#' (\code{is_marker}) in the column meta-data, indicating which columns are marker
#' columns; otherwise all columns are transformed.
#' 
#' 
#' @param d_se Input data. Assumed to be in the form of a 
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}}, prepared with the function 
#'   \code{\link{prepareData}}. Column meta-data is assumed to contain a vector of logical
#'   entries (\code{is_marker}) indicating marker columns.
#' 
#' @param cofactor Cofactor for 'arcsinh' transform. Default = 5, which is appropriate 
#'   for mass cytometry (CyTOF) data. For fluorescence flow cytometry, we recommend 
#'   cofactor = 150 instead.
#' 
#' 
#' @return d_se Data with transform applied to protein marker columns.
#' 
#' 
#' @importFrom SummarizedExperiment assays colData 'assays<-'
#' 
#' @export
#'
#' @examples
#' library(flowCore)
#' library(limma)
#' 
#' # filenames
#' files <- list.files(system.file("extdata", package = "diffcyt"), 
#'                     pattern = "\\.fcs$", full.names = TRUE)
#' files_BCRXL <- files[grep("BCRXL", files)]
#' files_ref <- files[grep("ref", files)]
#' 
#' # load data
#' files_load <- c(files_BCRXL, files_ref)
#' d_input <- lapply(files_load, read.FCS, transformation = FALSE, truncate_max_range = FALSE)
#' 
#' # sample IDs and group IDs
#' sample_IDs <- gsub("\\.fcs$", "", basename(files_load))
#' sample_IDs
#' group_IDs <- gsub("^patient[0-9]_", "", sample_IDs)
#' group_IDs
#' 
#' # indices of all marker columns, lineage markers, and functional markers
#' # (see Table 1 in Bruggner et al. 2014)
#' marker_cols <- c(3:4, 7:9, 11:19, 21:22, 24:26, 28:31, 33)
#' lineage_cols <- c(3:4, 9, 11,12,14, 21, 29, 31, 33)
#' functional_cols <- setdiff(marker_cols, lineage_cols)
#' 
#' # prepare data
#' d_se <- prepareData(d_input, sample_IDs, group_IDs, marker_cols, lineage_cols, functional_cols)
#' 
#' # transform data
#' d_se <- transformData(d_se, cofactor = 5)
transformData <- function(d_se, cofactor = 5) {
  
  is_marker <- colData(d_se)$is_marker
  
  if (is.null(marker_cols)) marker_cols <- 1:ncol(d_se[[1]])
  
  # extract expression data
  d_ex <- assays(d_se)[[1]]
  
  # transform marker columns
  d_ex[, is_marker] <- asinh(d_ex[, is_marker] / cofactor)
  
  assays(d_se)[[1]] <- d_ex
  
  d_se
}


