##########################################################################################
# Modified version of function 'tperm.fd()' from package 'fda' (Ramsay et al. 2014)
# 
# - Allow weighted permutation t-tests by accepting sample-specific weights in arguments 
# 'weights1' and 'weights2'
# - Use functions from 'matrixStats' package ('rowWeightedMeans', 'rowWeightedVars', 
# 'colMaxs') for faster runtime
# - Not calculating quantiles, point-wise p-values and q-values, plots; returning p-values
# only (for faster runtime)
# - Allow paired tests
# - Note: alternative version 'tperm.fd_fast.R' gives faster runtime when weights are not 
# required
#
# Modified by Lukas Weber and Mark Robinson, March 2017
##########################################################################################


#' @importFrom matrixStats rowWeightedMeans rowWeightedVars colMaxs
#' @importFrom fda is.fd eval.fd
#' 
.tperm.fd_wtd_fast <- function(x1fd, x2fd, weights1 = NULL, weights2 = NULL, 
                               paired = FALSE, nperm = 1000, argvals = NULL)
{
    if (is.null(weights1) | is.null(weights2)) {
        stop("use alternative version 'tperm.fd_fast.R' when weights are not required (faster runtime)")
    }
    if (!is.fd(x1fd) | !is.fd(x2fd)) {
        stop("x1fd and x2fd must both be functional data objects")
    }
    rangeobs = x1fd$basis$range
    rangehat = x2fd$basis$range
    if (!prod(rangeobs == rangehat)) {
        stop("x1fd and x2fd do not have the same range")
    }
    if (is.null(argvals)) {
        argvals = seq(rangeobs[1], rangeobs[2], length.out = 101)
    }
    x1mat = eval.fd(argvals, x1fd)
    x2mat = eval.fd(argvals, x2fd)
    n1 = ncol(x1mat)
    n2 = ncol(x2mat)
    weights = c(weights1, weights2)
    if (!((length(weights1) == n1) & (length(weights2) == n2))) {
        stop("lengths of weights vectors do not match number of samples")
    }
    Xmat = cbind(x1mat, x2mat)
    Tnullvals = matrix(0, length(argvals), nperm)
    if (!paired) {
        # unpaired (independent samples) permutation t-tests
        for (i in 1:nperm) {
            # keep track of column indices so each weight can be associated with the correct sample
            # weighted unequal-variance two-sample t-tests: 
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
        if (!(n1 == n2)) {
            stop("'n1' and 'n2' must be equal for paired permutation tests")
        }
        Xmat <- x2mat - x1mat
        weights <- weights1 + weights2  # paired tests: use combined number of cells as weights
        for (i in 1:nperm) {
            # paired permutation t-tests
            # https://stats.stackexchange.com/questions/64212/randomisation-permutation-test-for-paired-vectors-in-r
            # i.e. randomly flip signs (equivalent to permuting group labels within each pair)
            signs_i <- sample(c(-1, 1), n1, replace = TRUE)
            tXmat <- t(t(Xmat) * signs_i)  # transpose since multiplication wraps by column
            tmean <- rowWeightedMeans(tXmat, w = weights)
            tvar <- rowWeightedVars(tXmat, w = weights)
            Tnullvals[, i] <- abs(tmean) / sqrt(tvar/n1)  # note: dividing vars by n1 here
        }
        Tnull <- colMaxs(Tnullvals)
        means <- rowWeightedMeans(Xmat, w = weights)
        vars <- rowWeightedVars(Xmat, w = weights)
        Tvals <- abs(means) / sqrt(vars/n1)  # note: dividing vars by n1 here
        Tobs <- max(Tvals)
        pval <- mean(Tobs < Tnull)
    }
    return(list(pval = pval))
}


