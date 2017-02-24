#' Check input data format
#' 
#' Check input data is in the format required by 'diffcyt' functions
#' 
#' The functions in the 'diffcyt' package assume that input data is provided as as 
#' \code{\link[flowCore]{flowSet}} object (a set of \code{\link[flowCore]{flowFrame}} 
#' object) from the \code{\link[flowCore]{flowCore}} package.
#' 
#' This function checks whether input data is in the correct format, and returns an error
#' message if not.
#' 
#' 
#' @param d_input Input data.
#' 
#' 
#' @return d_input Input data, which has been confirmed as being in the correct format 
#'   (\code{\link[flowCore]{flowSet}} from the \code{\link[flowcore]{flowCore}} package).
#' 
#' 
#' @export
#'
#' @examples
#' # need to create a small example data set for examples
checkInputData <- function(d_input) {
  
  if (!is(d_input, "flowSet")) {
    stop(cat("Input data must be provided as a 'flowSet' object (a set of 'flowFrame'", 
             "objects) from the 'flowCore' package. See '?checkInputData' for an", 
             "example showing how to create a 'flowSet' from raw input data."))
  } else {
    print(TRUE)
    invisible(d_input)
  }

}


