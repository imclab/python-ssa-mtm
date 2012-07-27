# -*- coding: utf-8 -*-
#
# mssa.py
#
# purpose:  Follow the pdf "A beginner's guide to SSA".
# author:   Filipe P. A. Fernandes
# e-mail:   ocefpaf@gmail
# web:      http://ocefpaf.tiddlyspot.com/
# created:  12-Jul-2012
# modified: Fri 27 Jul 2012 05:32:36 PM BRT
#
# obs: 2. Multivariate Singular Spectral Analysis
#

from __future__ import division

import numpy as np
import scipy.stats.stats as stats
import matplotlib.pyplot as plt

#from scipy.linalg import toeplitz
#from oceans.ff_tools import lagcorr

# Some constants.
M, N = 4, 20
t = np.arange(N)


def gen_series(default=True):
    r"""The series is sampled from a sinus function with period 3.3, with
    Gaussian noise added, and has mean=0 and SD=1."""

    # Since I cannot reproduce the random seed for SciLab.
    # I just grabbed the values from the pdf.
    if default:
        X0 = [1.0135518, -0.7113242, -0.3906069, 1.565203, 0.0439317,
              -1.1656093, 1.0701692, 1.0825379, -1.2239744, -0.0321446,
              1.1815997, -1.4969448, -0.7455299, 1.0973884, -0.2188716,
              -1.0719573, 0.9922009, 0.4374216, -1.6880219, 0.2609807]

        X1 = [1.2441205, -0.7790017, 0.7752333, 1.4445798, -0.4838554,
              -0.8627650, 0.6981061, -0.7191011, -2.0409115, 0.3918897,
              1.1679199, -0.3742759, -0.1891236, 1.8180156, -0.4703280,
              -1.0197255, 0.9053288, -0.0818220, -1.3862948, -0.0379891]
    else:
        # Sine function with period length 3.3
        X0 = np.sin(2.0 * np.pi * t / 3.3)
        X1 = X0 ** 3
        # Original take a float, but numpy take only integers.
        np.random.seed(123456789)  # 0.123456789

        # Gaussian noise.
        noise = 0.2 * np.random.randn(N, 2)
        X0, X1 = X0 + noise[:, 0], X1 + noise[:, 0]
        # Remove mean value.
        X0, X1 = X0 - X0.mean(), X1 - X1.mean()
        # Normalize to std=1.
        X0, X1 = X0 / X0.std(), X1 / X1.std()

    return X0, X1


# Original series.
X0, X1 = gen_series()
fig, (ax0, ax1) = plt.subplots(nrows=2, ncols=1, figsize=(8, 6))
ax0.plot(X0, '-r.')
ax0.set_title('Time series X1(t) vs t')
ax1.plot(X1, '-r.')
ax1.set_title('Time series X1(t) vs t')
plt.show()


"""An essential and necessary step for MSSA is to normalize both time series.
That means to remove the mean value and to divide it by the standard deviation
(for each series separately)."""

X0_zs, X1_zs = stats.zscore(X0), stats.zscore(X1)


# Using shifted time series.
def shift(arr, n, order='forward'):
    if isinstance(arr, np.ndarray):
        arr = arr.tolist()
    if order == 'forward':
        shifted = arr[n:] + [0] * n
    elif order == 'reversed':
        shifted = [0] * n + arr[:-n]
    else:
        print("Order %s not recognized.  Try forward or reversed" % order)

    return shifted

# Embedded Time Series.
Y0 = np.c_[X0_zs, shift(X0_zs, 1), shift(X0_zs, 2), shift(X0_zs, 3)]
Y1 = np.c_[X1_zs, shift(X1_zs, 1), shift(X1_zs, 2), shift(X1_zs, 3)]
Y = np.c_[Y0, Y1]

# Covariance matrix.
C = np.dot(Y.T, Y) / N

# Computing the eigenvalues (lambda) eigenvectors (rho) or C:
lamb, rho = np.linalg.eig(C)
idx = lamb.argsort()[::-1]
rho = rho[:, idx]
lamb = lamb[idx]

rho[:, 0] = -rho[:, 0]  # Not sure why the pdf has this negative?

# First pair of eigenvectors.
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 6))
ax.set_title(r'Eigenvectors "$\rho$"')
ax.plot(rho[:, 0], 'b.-', label='1')
ax.plot(rho[:, 1], 'g.-', label='2')
ax.legend(numpoints=1)
plt.show()

# Principal components.
PC = np.dot(Y, rho)

"""By writing down element-by-element the results of this matrix product (try
this!) one discovers that each element of the PC matrix is the sum of a linear
combination of M values of the first time series (weighted by that part of the
EOF that corresponds to the first time series) and a linear combination of the
M values of the second time series (again weighted by the corresponding part of
the EOF).  This means that each PC contains characteristics of both time
series.  Unlike the EOFs or the matrices Y and C, we can no longer identify a
part that corresponds to each separate time series."""

fig, (ax0, ax1, ax2, ax3) = plt.subplots(nrows=4, ncols=1, figsize=(8, 6))
ax0.set_title(r'Principal components PC vs t')
ax0.plot(PC[:, 0], 'b.-')
ax1.plot(PC[:, 1], 'b.-')
ax2.plot(PC[:, 2], 'b.-')
ax3.plot(PC[:, 3], 'b.-')
plt.show()

# Reconstruction of the time series.
RC0, RC1 = np.zeros((N, M)), np.zeros((N, M))

for m in np.arange(M):
    Z = np.zeros((N, M))  # Time-delayed embedding of PC[:, m].
    for m2 in np.arange(M):
        Z[m2 - N:, m2] = PC[:N - m2, m]

    # Determine RC as a scalar product.
    RC0[:, m] = np.dot(Z, rho[:M, m] / M)
    RC1[:, m] = np.dot(Z, rho[M:2 * M, m] / M)

fig, axs = plt.subplots(nrows=4, ncols=2, figsize=(8, 6))
fig.suptitle(r'Reconstruction components RC vs t')
axs[0, 0].plot(RC0[:, 0], 'r.-')
axs[1, 0].plot(RC0[:, 1], 'r.-')
axs[2, 0].plot(RC0[:, 2], 'r.-')
axs[3, 0].plot(RC0[:, 3], 'r.-')

axs[0, 1].plot(RC1[:, 0], 'r.-')
axs[1, 1].plot(RC1[:, 1], 'r.-')
axs[2, 1].plot(RC1[:, 2], 'r.-')
axs[3, 1].plot(RC1[:, 3], 'r.-')
plt.show()

"""The first both RC1 and RC2 describe the oscillations, where the RC3 and RC4
describe a trend (which may be introduced due to the random number generator
and the very short time series).  When we summarize all eight RCs we
reconstruct the whole time series."""

fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8, 6))
fig.suptitle(r'Original time series and reconstructions vs t')
axs[0, 0].plot(X0, 'b-', label='X0 Original')
axs[0, 0].plot(RC0.sum(axis=1), 'r.-', label='RC0-7')
axs[0, 0].legend(numpoints=1)

axs[1, 0].plot(X1, 'b-', label='X1 Original')
axs[1, 0].plot(RC1.sum(axis=1), 'r.-', label='RC0-7')
axs[1, 0].legend(numpoints=1)

axs[0, 1].plot(X0, 'b-', label='X0 Original')
axs[0, 1].plot(RC0[:, 0] + RC0[:, 1], 'r.-', label='RC0-1')
axs[0, 1].legend(numpoints=1)

axs[1, 1].plot(X1, 'b-', label='X1 Original')
axs[1, 1].plot(RC1[:, 0] + RC1[:, 1], 'r.-', label='RC0-1')
axs[1, 1].legend(numpoints=1)

plt.show()

"""Advantages of MSSA with respect to SSA:
The MSSA allows in the same way as SSA to decompose the time series into its
spectral components.  Like in single-variate SSA, we are thus able to identify
trends and oscillating pairs.  But in contrast to SSA, the MSSA also takes
cross-correlations into account, where MSSA is a combination of SSA and
principal component analysis (PCA).  The individual RCs of the different time
series are connected; they represent the same spectral part.  We are hence
able to identify oscillatory components (e.g. limit cycles) that are intrinsic
to all time series."""
