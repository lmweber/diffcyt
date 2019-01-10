#' Alias for function 'topTable' (deprecated)
#' 
#' Alias for function 'topTable' (deprecated)
#' 
#' The function \code{topClusters} has been renamed to \code{\link{topTable}}, to more
#' accurately reflect the structure of the results (results are returned for either
#' clusters or cluster-marker combinations, depending on the type of differential tests
#' performed).
#' 
#' This alias is provided for backward compatibility. The new function name
#' \code{\link{topTable}} should be used whenever possible.
#' 
#' See \code{\link{topTable}} for details.
#' 
#' @param ... See arguments for function \code{\link{topTable}}
#' 
#' @export
#' 
topClusters <- function(...) {
  .Deprecated("topTable")
  topTable(...)
}


