# -*- coding: utf-8 -*-
#
# ssa.py
#
# purpose:  Reproduce Julien's hepta_ssa.m in python
# author:   Filipe P. A. Fernandes
# e-mail:   ocefpaf@gmail
# web:      http://ocefpaf.tiddlyspot.com/
# created:  23-Jul-2012
# modified: Fri 27 Jul 2012 05:42:54 PM BRT
#
# obs:
#

from __future__ import division

import numpy as np
from scipy import optimize
from scipy.linalg import toeplitz
from matplotlib.mlab import prctile
from scipy.stats import nanstd, nanmean

eps = np.spacing(1)


def xcorr(x, y=None, matlab_compat='unbiased', maxlags=None):
    r"""Reproduce Matlab's xcorr behavior.
    Adapted from matplotlib.pyplot import xcorr."""

    Nx = len(x)
    if not y:
        y = x.copy()
    else:
        if Nx != len(y):
            raise ValueError('x and y must be equal length.')

    c = np.correlate(x, y, mode='full')

    if maxlags is None:
        maxlags = Nx - 1

    if maxlags >= Nx or maxlags < 1:
        raise ValueError('maglags must be None or positive < %d.' % Nx)

    lags = np.arange(-maxlags, maxlags + 1)
    c = c[Nx - 1 - maxlags:Nx + maxlags]

    # TODO: Add all the other options.
    if matlab_compat == 'unbiased':
        c /= (Nx - np.abs(lags))

    return c, lags


def standardize(X, axis=None):
    r"""Subtracts the data mean and divide by its standard deviation
    at the specified axis.  Accepts NaNs."""

    # NOTE: There is an alternative in scipy.stats.mstats.zscore
    mu = nanmean(X, axis=axis)
    sigma = nanstd(X, axis=axis)

    Xr = (X - mu) / sigma

    return Xr, mu, sigma


def eigd(A):
    r"""Eigenvalues in decreasing order, and corresponding eigenvectors.

    V, d = eigd(A) determines the positive eigenvalues `d` and corresponding
    eigenvectors `V` of a matrix `A`.  The eigenvalues are returned as the
    vector `d`.  The eigenvectors are returned as the columns of the matrix
    `V`.  The function `eigd` calls the function eig, which contrary to `eigs`
    computes all eigenvalues and eigenvectors.

    See also: `peigs`.
    """

    m, n = A.shape

    # Eigendecomposition of A.
    d, V = np.linalg.eig(A)

    # Ensure that eigenvalues are monotonically decreasing.
    idx = d.argsort()[::-1]
    V = V[:, idx]
    d = d[idx]

    return V, d


def pca_truncation_criteria(d, p, r, N, b=1):
    r"""Computes criteria for truncating principal component analyses.

    wk85, ne08 = pca_truncation_criteria(d, p, r, N, b) returns the truncation
    choice criteria Mean Description Length, Akaike's Information Criterion,
    and the bias-corrected AIC, given as inputs a vector of eigenvalues `d`,
    the problem dimension `p`, possible truncation levels `r`, and sample size
    `N`.  The optional argument `b` indicates whether the problem is real
    valued (b=1) or complex-valued (b=2); the default is b=1.

    The implementation of the truncation choice criteria follows Wax & Kailath,
    "Detection of Signals by Information Theoretic Criteria," IEEE Trans.
    Acoust. Speech Signal Proc., 33 (1985), 387--392. It uses Wax & Kailath's
    (WK85) or Nadakuditi & Edelman (NE08) log-likelihood function and  with
    that computes various information-theoretic criteria to select the
    truncation level:

    MDL: Schwartz (1978) and Rissanen (1978)
    AIC: Akaike (1973, 1974) [times 1/2 to obtain likelihood measure]
    AICC: Hurvich & Tsai (1989) [times 1/2]

    If the problem is rank-deficient (number of nonzero eigenvalues d < p), it
    uses an ad hoc restricted log-likelihood.  In those cases, NE08 may be
    the preferred criterion (as it is designed to handle such cases).

    Outputs, structures wk85 and ne08, with vectors mdl, aic and aicc as
    fields.

    Tapio Schneider, 2/18/2012
    J.E.G. added BIC, AICC to NE08, which was a modified AIC.
    """
    d, p, r, N, b = map(np.asanyarray, (d, p, r, N, b))
    # Sort eigenvalues (in case they are not).
    idx = d.argsort()[::-1]
    d = d[idx]

    # Number of numerically nonzero eigenvalues,
    d_min = eps * d.max() * p
    nd = np.sum(d > d_min)

    rmax = (np.r_[r.max(), len(d) - 1]).min()
    r = r[0:rmax]

    if r.max() > len(d) - 1:
        raise ValueError('Truncation level too large.')

    # Assemble log-likelihood (Wax & Kailath 1985) and Nadakuditi & Edelman
    # (2008) statistic for various truncations r.
    llh = np.zeros_like(r, dtype=np.float)
    tr = np.zeros_like(r, dtype=np.float)
    for j in range(0, len(r)):
        k = r[j]
    # Wax & Kailath log-likelihood (ad hoc restricted to nonsingular subspace
    # if nd < p, i.e., if sample covariance matrix is singular).
    llh[j] = N * (nd - k) * np.log(np.prod(d[k:nd]) ** (1.0 / (p - k)) /
                                   np.sum(d[k:nd]) * (p - k))

    # Nadakuditi & Edelman statistic.
    tr[j] = p * (((p - k) * np.sum(d[k:nd] ** 2) /
                 np.sum(d[k:nd]) ** 2.0 - (1.0 + p / N)) -
                 (2.0 / b - 1.0) * p / N)

    """Number of free parameters on which log-likelihood depends (number of
    eigenvector and eigenvalue components, minus the number of constraints for
    orthogonality and unit norm of eigenvectors).

    This number is obtained as follows:

    We have k + 1 + b * r * p parameters for eigenvalues, eigenvectors, and
    noise variance. From that subtract the number of constraints:
    -- normalization: b[(r - 1) + (kr1) + ... + 1] = br / 2(r - 1) constraints
    -- orthogonality: b * r constraints
    Total: r * (bp - (b / 2) r - b/2 + 1) + 1 free parameters."""
    peff = r * (b * p - b / 2.0 * r - b / 2 + 1) + 1

    # Assemble truncation choice criteria
    wk85 = dict()
    mdl = -llh + peff * np.log(N) / 2.0
    # For AIC(c), use 1/2 the traditional measure to obtain a likelihood
    # measure.
    aic = -llh + peff
    aicc = aic + peff * (peff + 1.0) / (N - peff - 1)

    wk85.update(dict(mdl=mdl, aic=aic, aicc=aicc))

    ne08 = dict()
    pp = r + 1  # Number of paramters.
    aic = b / 8.0 * (N / p) ** 2.0 * tr ** 2 + pp
    aicc = aic + peff * (pp + 1) / (N - pp - 1)
    mdl = b / 8.0 * (N / p) ** 2 * tr ** 2 + pp * np.log(N) / 2.0

    ne08.update(dict(mdl=mdl, aic=aic, aicc=aicc))

    return wk85, ne08


def ar1(x):
    r"""Allen and Smith AR(1) model estimation.
    Syntax: g, a, mu2 = ar1(x)

    Input:  x - time series (univariate).

    Output: g - estimate of the lag-one autocorrelation.
        a - estimate of the noise variance.
        mu2 - estimated square on the mean.

    AR1 uses the algorithm described by Allen and Smith 1995, except that
    Matlab's 'fzero' is used rather than Newton-Raphson.

    Fzero in general can be rather picky - although I haven't had any problem
    with its implementation here, I recommend occasionally checking the output
    against the simple estimators in AR1NV.

    Alternative AR(1) estimatators: ar1cov, ar1nv, arburg, aryule

    Written by Eric Breitenberger.      Version 1/21/96
    Please send comments and suggestions to eric@gi.alaska.edu
    Updated,optimized&stabilized by Aslak Grinsted 2003-2005.
    """

    N = len(x)
    m = x.mean()
    x = x - m

    # Lag zero and one covariance estimates:
    c0 = np.dot(x, x / N)
    c1 = np.dot(x[0:N - 1], x[1:N]) / (N - 1)

    g0 = c1 / c0  # Initial estimate for gamma.

    # Optimize gammest.
    def gammest(gin):
        r"""Used by AR1 to compute a function for minimization by fzero.

        Written by Eric Breitenberger.  Version 1/21/96
        Please send comments and suggestions to eric@gi.alaska.edu
        """

        gk = np.arange(1, N)
        gk = gin ** gk
        mu2 = (1.0 / N) + (2.0 / N ** 2.0) * np.sum(np.arange(N - 1, 0, -1) *
                                                    gk)
        gout = (1.0 - g0) * mu2 + g0 - gin
        if gout > 1:
            gout = np.NaN

        return gout

    # Find g by getting zero of `gammest`.
    # There are tons of optimizition algorithms in SciPy.  I'm not sure which
    # compares better with fzero.
    g = optimize.newton(gammest, g0, tol=0.0001)

    gk = np.arange(1, N)
    gk = g ** gk
    mu2 = (1.0 / N) + (1.0 / N ** 2.0) * 2.0 * np.sum(np.arange(N - 1, 0, -1) *
                                                      gk)
    c0est = c0 / (1.0 - mu2)
    a = np.sqrt((1.0 - g ** 2) * c0est)

    return g, a, mu2


def ssa(X, M=None, K=0):
    r"""Performs Singular Spectrum Analysis on time series X with the method of
    Vautard and Ghil, Phys. D. 1989.

    Parameters
    ----------
    X : 1D array
        Vector of evenly spaced observations.
    M : int
        Window length.  Default value is M = len(X) / 10
    K : int
        Number of EOFs used for reconstruction (AICC choice by default k=0).
        if K = 0, corrected Akaike Information Criterion (AICC) is used
        if K = 'mcssa', the Monte Carlo spectral significance estimation of
        Allen & Smith (J Clim, 1996) is used.

    Returns
    -------
    spec : array_like
           Eigenvalue spectrum, in % variance.
    eig_vec : array_like
              Eigenvector matrix ("temporal EOFs").
    PC : array_like
         Matrix of principal components.
    RC : array_like
         Matrix of RCs (N*M, K) (only if K > 0).
    RCp : array_like
          Reconstructed time-series, involving only the modes retained, and
          rescaled to original mean and variance.

    Examples
    --------
    spec, eig_vec, PC, RC, RCp = ssa(X,[M, K])

    Notes
    -----
    Orignal file hepta_ssa.m from Hepta Technologies, 2004 writing in MatlabTM.
    last updated 03/14/2012 to include automated choice for K (AICC).

    Julien Emile-Geay, Lamont Doherty Earth Observatory.
    Dec 2004 last updated 03/14/2012
    """

    X = np.atleast_1d(X)
    if X.ndim > 1:
        raise ValueError("Input vector `X` has more than 1 dimension.")

    N = len(X)

    # Center the series.
    Xr, mu, sigma = standardize(X)  # NOTE: Original calls standardize.m.

    # Set default value for M.
    if not M:
        M = N // 10
    if K == 'mcssa':
        mcssa = True
        MC = 1000
    else:
        mcssa = False
        signif = np.arange(0, K)  # FIXME: 0, K

    Np = N - M + 1

    gam, lags = xcorr(Xr, maxlags=M - 1, matlab_compat='unbiased')

    # Fill in Covariance matrix.  Take positive half of auto-correlation
    # diagram, hence M to 2M - 1.
    C = toeplitz(gam[M - 1:2 * M])

    # Solve eigenvalue problem.
    eig_vec, eig_val = eigd(C)  # FIXME: Matlab eig_vec have reversed signs.
    spec = eig_val / np.sum(eig_val)

    # Determine significant eigenvalues.
    if mcssa:
        # NOTE: Got this at from: http://www.gps.caltech.edu/~tapio/arfit/
        # But this is commented out in the original code.
        #w, A, C, SBC, FPE, th = arfit(Xr, 1, 1)  # fit AR(1) model.
        # NOTE: The original code uses ar1.m.
        # What is the difference between ar1.m and arfit.m?
        a, var, _ = ar1(Xr)
        s = np.sqrt(var)
        noise = np.zeros(N, MC)
        noise[0, :] = np.tile(Xr[0], np.r_[1, MC])
        for jt in range(1, N):
            noise[jt, :] = a * noise[jt - 1, :] + s * np.random.randn(1, MC)

        noise, _, _ = standardize(noise)
        Lambda_R = np.zeros_like(MC)  # FIXME: Don't know the right shape yet.
        for m in range(0, MC):
            Gn, ln = xcorr(noise[:, m], M - 1, 'unbiased')
            Cn = toeplitz(Gn[M: 2 * M - 1])
            # Noise "eigenvalues".
            tmp = np.dot(eig_vec, Cn)
            Lambda_R[:, m] = np.diag(np.dot(tmp, eig_vec))

        q95 = prctile(Lambda_R, 100 * 0.95)  # FIXME
        # Index of modes rising above the background.
        signif = np.where(eig_val > q95)
        print('MCSSA modes retained: %s' % signif)

        fix, ax = plt.subplots()
        ax.set_title('MCSSA')
        v = np.arange[1, M + 1]
        ligr = [0.7000, 0.7000, 0.7000]
        lmin = Lambda_R.min(axis=1)
        lmax = Lambda_R.max(axis=1)
        ax.fill(v, lmin, lmax, ligr, ligr, 0, 0.3)
        ax.plot(v, eig_val, 'kx', linewidth=2.0)
        ax.plot(v, q95, 'r-', linewidth=2.0)
    elif K == 0:
        trunc = range(0, len(spec))
        # The pca_truncation_criteria.m original call:
        # [MDL, NE08, AIC, AICC] =
        # pca_truncation_criteria(eig_val, 1, trunc, N, 1)
        WK85, NE08 = pca_truncation_criteria(eig_val, 1, trunc, N, 1)
        imin = (np.real(NE08['aicc'])).argmin()
        K = trunc[imin]
        print('AICC truncation choice, K = %s' % K)
        signif = np.arange(0, K)

    # Compute PCs.
    decal = np.zeros((Np, M))

    for t in range(0, N - M + 1):
        decal[t, :] = Xr[t:M + t]

    # The columns of this matrix are Ak(t), k=1 to M.
    PC = np.dot(decal, eig_vec)

    # Compute reconstructed timeseries if K > 0.
    if signif:
        RC = np.zeros((N, len(signif)))
        # First M terms.
        for t in range(0, M - 1):
            Av = np.flipud(PC[0:t, signif])
            eig_vec_red = eig_vec[0:t, signif]
            RC[t, :] = 1.0 / t * np.sum(Av * eig_vec_red, axis=0)

        # Middle of timeseries.
        for t in range(M, Np + 1):
            Av = np.flipud(PC[t - M + 1:t, signif])
            eig_vec_red = eig_vec[0:M, signif]
            RC[t, :] = 1 / M * np.sum(Av * eig_vec_red, axis=0)

        # Last M terms.
        for t in range(Np + 1, N + 1):
            Av = np.flipud(PC[t - M + 1:Np, signif])
            eig_vec_red = eig_vec[t - N + M:M, signif]
            RC[t, :] = 1.0 / (N - t + 1) * np.sum(Av * eig_vec_red, axis=0)

        # Sum and restore the mean and variance.
        RCp = sigma * np.sum(RC, axis=1) + mu
    else:
        RC, RCp = None, None

    return spec, eig_vec, PC, RC, RCp, signif

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    X = [1.0135518, -0.7113242, -0.3906069, 1.565203, 0.0439317, -1.1656093,
         1.0701692, 1.0825379, -1.2239744, -0.0321446, 1.1815997, -1.4969448,
         -0.7455299, 1.0973884, -0.2188716, -1.0719573, 0.9922009, 0.4374216,
         -1.6880219, 0.2609807]

    spec, eig_vec, PC, RC, RCp, signif = ssa(X)

    # Plotting the PCs.
    fig, (ax0, ax1) = plt.subplots(nrows=2, ncols=1, figsize=(8, 6))
    ax0.set_title(r'Principal components PC vs T')
    ax0.plot(PC[:, 0], '.-')
    ax1.plot(PC[:, 1], '.-')
