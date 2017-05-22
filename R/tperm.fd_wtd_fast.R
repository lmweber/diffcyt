##########################################################################################
# Modified version of function 'tperm.fd()' from package 'fda' (Ramsay et al. 2014)
# 
# Modifications:
# - Allow sample-specific weights in arguments 'weights1' and 'weights2' (i.e. 'weights1' 
# and 'weights2' contain weights by sample in functional data objects 'x1fd' and 'x2fd'), 
# for weighted permutation t-tests.
# - Use 'rowWeightedMeans', 'rowWeightedVars', and 'colMaxs' from the 'matrixStats'
# package for faster runtime.
# - Not calculating quantiles, point-wise p-values / q-values, plots; returning p-values
# only (for faster runtime).
# - Allow paired tests.
# 
# Note: Alternative version 'tperm.fd_fast.R' gives faster runtime when weights are not 
# required.
#
# Modified by Lukas Weber and Mark Robinson, March 2017
##########################################################################################


#' @importFrom matrixStats rowWeightedMeans rowWeightedVars colMaxs
#' @importFrom fda is.fd eval.fd
#' @importFrom stats quantile
#' @importFrom graphics plot lines abline legend
#' 
.tperm.fd_wtd_fast <- function (x1fd, x2fd, weights1 = NULL, weights2 = NULL, paired = FALSE, 
                                nperm = 1000, q = 0.05, argvals = NULL, plotres = TRUE, ...)
  # note: leaving unused arguments (q, plotres, ...) for consistency with function signature from 'fda' package
{
    if (!is.fd(x1fd) | !is.fd(x2fd)) {
        stop("x1fd and x2fd must both be functional data objects")
    }
    rangeobs = x1fd$basis$range
    rangehat = x2fd$basis$range
    if (!prod(rangeobs == rangehat)) {
        stop("x1fd and x2fd do not have the same range.")
    }
    if (is.null(argvals)) {
        argvals = seq(rangeobs[1], rangeobs[2], length.out = 101)
    }
    x1mat = eval.fd(argvals, x1fd)
    x2mat = eval.fd(argvals, x2fd)
    n1 = ncol(x1mat)
    n2 = ncol(x2mat)
    weights = c(weights1, weights2)
    if (!is.null(weights)) {
        if (!((length(weights1) == n1) & (length(weights2) == n2))) {
            stop("lengths of weights vectors do not match number of samples")
        }
    }
    if (is.null(weights)) {
        warning(paste0("Use alternative version 'tperm.fd_fast.R' for faster runtime ", 
                       "when weights are not required."))
        weights1 = rep(1, n1)
        weights2 = rep(1, n2)
        weights = c(weights1, weights2)
    }
    Xmat = cbind(x1mat, x2mat)
    Tnullvals = matrix(0, length(argvals), nperm)
    if (!paired) {
        # unpaired (independent samples) permutation t-tests
        for (i in 1:nperm) {
            # keep track of column indices so each weight can be associated with the correct sample
            # formulas for weighted unequal-variance two-sample t-tests: 
            # http://influentialpoints.com/Training/two-sample_t-test-principles-properties-assumptions.htm
            perm_i = sample(n1 + n2)
            weights_i = weights[perm_i]
            tXmat = Xmat[, perm_i]
            tmean1 = rowWeightedMeans(tXmat, w = weights_i, cols = 1:n1)  # note: 'cols' argument subsets both x and w
            tmean2 = rowWeightedMeans(tXmat, w = weights_i, cols = n1 + 1:n2)
            tvar1 = rowWeightedVars(tXmat, w = weights_i, cols = 1:n1)
            tvar2 = rowWeightedVars(tXmat, w = weights_i, cols = n1 + 1:n2)
            Tnullvals[, i] = abs(tmean1 - tmean2) / sqrt(tvar1/n1 + tvar2/n2)  # note: dividing vars by n1 and n2 here
        }
        Tnull = colMaxs(Tnullvals)
        mean1 = rowWeightedMeans(Xmat, w = weights, cols = 1:n1)
        mean2 = rowWeightedMeans(Xmat, w = weights, cols = n1 + 1:n2)
        var1 = rowWeightedVars(Xmat, w = weights, cols = 1:n1)
        var2 = rowWeightedVars(Xmat, w = weights, cols = n1 + 1:n2)
        Tvals = abs(mean1 - mean2) / sqrt(var1/n1 + var2/n2)  # note: dividing vars by n1 and n2 here
        Tobs = max(Tvals)
        pval = mean(Tobs < Tnull)
    } else {
        # paired permutation t-tests
        Xmat <- x2mat - x1mat
        weights <- weights1 + weights2  # paired tests: use combined number of cells as weights
        for (i in 1:nperm) {
            # paired permutation t-tests
            # https://stats.stackexchange.com/questions/64212/randomisation-permutation-test-for-paired-vectors-in-r
            # (i.e. randomly flip sign of each paired difference)
            signs_i <- sample(c(-1, 1), length(argvals), replace = TRUE)
            tXmat <- t(t(Xmat) * signs_i)  # transpose since multiplication wraps by column
            tmean <- rowWeightedMeans(tXmat, w = weights)
            tvar <- rowWeightedVars(tXmat, w = weights)
            stopifnot(n1 == n2)
            Tnullvals[, i] <- abs(tmean) / sqrt(tvar/n1)
        }
        Tnull <- colMaxs(Tnullvals)
        means <- rowWeightedMeans(Xmat, w = weights)
        vars <- rowWeightedVars(Xmat, w = weights)
        stopifnot(n1 == n2)
        Tvals <- abs(means) / sqrt(vars/n1)
        Tobs <- max(Tvals)
        pval <- mean(Tobs < Tnull)
    }
    return(list(pval = pval))
}


