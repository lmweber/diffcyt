linear_sum <- function(covariates_list,betas, random_covariates_list = NULL){
  stopifnot(is.list(covariates_list))
  stopifnot(is.list(random_covariates_list) | is.null(random_covariates_list))
  stopifnot(length(betas)==length(covariates_list))
  # check for equal list lengths
  stopifnot(length(unique(purrr::map(covariates_list, ~ length(.x)))) == 1)
  stopifnot((length(unique(purrr::map(random_covariates_list, ~ length(.x)))) == 1) | is.null(random_covariates_list))
  lin_sum <- matrix(unlist(covariates_list),ncol = length(covariates_list)) %*% matrix(betas,ncol=1)
  if (!is.null(random_covariates_list)){
    lin_sum <- lin_sum +
      matrix(unlist(random_covariates_list), ncol = length(random_covariates_list)) %*%
      matrix(rep(1,length(random_covariates_list)),ncol=1)
  }
  return(lin_sum)
}

# if multi cluster are simulated each batch of one differential cluster plus
# multiple non-differential clusters add up to one which needs to change,
# so that all batches (all clusters) together add to one 
calculate_wanted_cluster_proportion <- function(real_sizes){
  cur_prop <- sum(real_sizes)
  stopifnot(cur_prop <= 1)
  wanted_prop <- real_sizes/cur_prop
  increase <- 1/cur_prop
  stopifnot(increase>=2)
  stopifnot(all.equal(sum(wanted_prop),1))
  return(wanted_prop)
}

is_valid_data_transform <- function(transform_fn){
  return(10 == length(transform_fn(seq_len(10))))
}


boxcox_transform <- function(x,formula,verbose=FALSE){
  bc <- MASS::boxcox(formula, plotit=FALSE)
  lambda <- bc$x[which.max(bc$y)]
  if (verbose) cat(paste0("boxcox-transform lambda = ",lambda,"\n"))
  if (lambda==0){
    return(log(x))
  } else {
    return((x^(lambda)-1)/lambda)  }
}

#' @importFrom stats rnorm
do_transformation <- function(x,transform_fn,verbose=FALSE){
  if (identical(transform_fn,"identity")) {
    return(x)
  } else if (identical(transform_fn,"boxcox")) {
    tmp <- rnorm(length(x))
    return(boxcox_transform(x,x~tmp, verbose))
  } else if (identical(transform_fn,"boxcox_positive")) {
    tmp <- rnorm(length(x))
    if(verbose) cat("\nTrVal:\t")
    x <- boxcox_transform(x,x~tmp, verbose)
    if (min(x)<=0){
      x <- x-min(x)*1.00001 # ensure non zero values for ppd
    }
    return(x)
  } else if (identical(transform_fn,"log_positive")) {
    x <- log(x)
    if (min(x)<=0){
      x <- x-min(x)*1.00001 # ensure non zero values for ppd
    }
    return(x)
  } else if (is(transform_fn, "character")){
    return(tryCatch(do.call(transform_fn, list(x=x)),
                   error = function(e) {
                     stop("'transform_fn' is no valid transformation",call. = FALSE)
                   }))
  } else if (is_valid_data_transform(transform_fn)){
    return(transform_fn(x))
  } else {
    stop("'transform_fn' is no valid transformation",call. = FALSE)
  }
}


