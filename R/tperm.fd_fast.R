##########################################################################################
# Modified version of function 'tperm.fd()' from package 'fda' (Ramsay et al. 2014)
# 
# - Use 'rowMeans' and functions from 'matrixStats' package ('rowVars', 'colMaxs') for 
# faster runtime
# - Not calculating quantiles, point-wise p-values and q-values, plots; returning p-values
# only (for faster runtime)
# - Allow paired tests
# - Note: alternative version 'tperm.fd_wtd_fast.R' also accepts sample-specific weights
# to allow weighted tests, but has slower runtime
# 
# Modified by Lukas Weber and Mark Robinson, March 2017
##########################################################################################


#' @importFrom matrixStats rowVars colMaxs
#' @importFrom fda is.fd eval.fd
#' 
.tperm.fd_fast <- function(x1fd, x2fd, paired = FALSE, nperm = 1000, argvals = NULL)
{
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
    if (!paired) {
      # unpaired (independent samples) permutation t-tests
        Xmat = cbind(x1mat, x2mat)
        Tnullvals = matrix(0, length(argvals), nperm)
        for (i in 1:nperm) {
            tXmat = Xmat[, sample(n1 + n2)]
            tmean1 = rowMeans(tXmat[, 1:n1])
            tmean2 = rowMeans(tXmat[, n1 + (1:n2)])
            tvar1 = rowVars(tXmat, cols = 1:n1)
            tvar2 = rowVars(tXmat, cols = n1 + 1:n2)
            Tnullvals[, i] = abs(tmean1 - tmean2)/sqrt(tvar1/n1 + tvar2/n2)  # note: dividing vars by n1 and n2 here
        }
        Tnull = colMaxs(Tnullvals)
        mean1 = rowMeans(Xmat[, 1:n1])
        mean2 = rowMeans(Xmat[, n1 + (1:n2)])
        var1 = rowVars(Xmat, cols = 1:n1)
        var2 = rowVars(Xmat, cols = n1 + 1:n2)
        Tvals = abs(mean1 - mean2)/sqrt(var1/n1 + var2/n2)  # note: dividing vars by n1 and n2 here
        Tobs = max(Tvals)
        pval = mean(Tobs < Tnull)
    } else {
        # paired permutation t-tests
        if (!(n1 == n2)) {
            stop("'n1' and 'n2' must be equal for paired permutation tests")
        }
        Xmat <- x2mat - x1mat
        for (i in 1:nperm) {
            # randomly flip signs (equivalent to permuting group labels within each pair)
            signs_i <- sample(c(-1, 1), n2, replace = TRUE)
            tXmat <- t(t(Xmat) * signs_i)  # transpose since multiplication wraps by column
            tmean <- rowMeans(tXmat)
            tvar <- rowVars(tXmat)
            Tnullvals[, i] <- abs(tmean) / sqrt(tvar/n1)  # note: dividing vars by n1 here
        }
        Tnull <- colMaxs(Tnullvals)
        means <- rowMeans(Xmat)
        vars <- rowVars(Xmat)
        Tvals <- abs(means) / sqrt(vars/n1)  # note: dividing vars by n1 here
        Tobs <- max(Tvals)
        pval <- mean(Tobs < Tnull)
    }
    return(list(pval = pval))
}

