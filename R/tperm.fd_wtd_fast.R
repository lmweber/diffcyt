##########################################################################################
# Modified version of function 'tperm.fd()' from package 'fda' (Ramsay et al. 2014)
# 
# Two main modifications:
# - Allow sample-specific weights in arguments 'weights1' and 'weights2' (i.e. 'weights1' 
# and 'weights2' contain weights by sample in functional data objects 'x1fd' and 'x2fd' 
# respectively), for weighted permutation t-tests.
# - Use 'rowWeightedMeans', 'rowWeightedVars', and 'colMaxs' from the 'matrixStats'
# package for faster runtime.
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
.tperm.fd_wtd_fast <- function (x1fd, x2fd, weights1 = NULL, weights2 = NULL, nperm = 200, 
                                q = 0.05, argvals = NULL, plotres = TRUE, ...)
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
    q = 1 - q
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
    for (i in 1:nperm) {
        # keep track of indices so each weight can be associated with the correct sample
        perm_i = sample(n1 + n2)
        weights_i = weights[perm_i]
        tXmat = Xmat[, perm_i]
        tmean1 = rowWeightedMeans(tXmat[, 1:n1], weights_i[1:n1])
        tmean2 = rowWeightedMeans(tXmat[, n1 + (1:n2)], weights_i[n1 + (1:n2)])
        tvar1 = rowWeightedVars(tXmat[, 1:n1], weights_i[1:n1]) / sum(weights_i[1:n1])
        tvar2 = rowWeightedVars(tXmat[, n1 + (1:n2)], weights_i[n1 + (1:n2)]) / sum(weights_i[n1 + (1:n2)])
        Tnullvals[, i] = abs(tmean1 - tmean2)/sqrt(tvar1 + tvar2)
    }
    Tnull = colMaxs(Tnullvals)
    mean1 = rowWeightedMeans(Xmat[, 1:n1], weights1)
    mean2 = rowWeightedMeans(Xmat[, n1 + (1:n2)], weights2)
    var1 = rowWeightedVars(Xmat[, 1:n1], weights1) / sum(weights1)
    var2 = rowWeightedVars(Xmat[, n1 + (1:n2)], weights2) / sum(weights2)
    Tvals = abs(mean1 - mean2)/sqrt(var1 + var2)
    Tobs = max(Tvals)
    pval = mean(Tobs < Tnull)
    qval = quantile(Tnull, q)
    pvals.pts = rowMeans(Tvals < Tnullvals)
    qvals.pts = apply(Tnullvals, 1, quantile, q)
    if (plotres) {
        if (is.null(names(x1fd$fdnames)) | is.null(names(x2fd$fdnames))) {
            xlab = "argvals"
        }
        else if (prod(names(x1fd$fdnames)[1] == names(x2fd$fdnames)[1])) {
            xlab = names(x1fd$fdnames)[1]
        }
        else {
            xlab = "argvals"
        }
        ylims = c(min(Tvals, qvals.pts), max(Tobs, qval))
        plot(argvals, Tvals, type = "l", col = 2, ylim = ylims, 
            lwd = 2, xlab = xlab, ylab = "t-statistic", ...)
        lines(argvals, qvals.pts, lty = 3, col = 4, lwd = 2)
        abline(h = qval, lty = 2, col = 4, lwd = 2)
        legendstr = c("Observed Statistic", paste("pointwise", 
            1 - q, "critical value"), paste("maximum", 1 - q, 
            "critical value"))
        legend(argvals[1], ylims[2], legend = legendstr, col = c(2, 
            4, 4), lty = c(1, 3, 2), lwd = c(2, 2, 2))
    }
    return(list(pval = pval, qval = qval, Tobs = Tobs, Tnull = Tnull, 
        Tvals = Tvals, Tnullvals = Tnullvals, qvals.pts = qvals.pts, 
        pvals.pts = pvals.pts, argvals = argvals))
}

