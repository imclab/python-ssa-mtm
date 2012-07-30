# -*- coding: utf-8 -*-
#
# ssa_tutorial.py
#
# purpose:  Follow the pdf "A beginner's guide to SSA".
# author:   Filipe P. A. Fernandes
# e-mail:   ocefpaf@gmail
# web:      http://ocefpaf.tiddlyspot.com/
# created:  29-Jun-2012
# modified: Fri 27 Jul 2012 05:43:09 PM BRT
#
# obs: The original pdf uses SciLab instead of python.
#

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

from scipy.linalg import toeplitz

from oceans.ff_tools import lagcorr


def gen_series(default=True):
    r"""The series is sampled from a sinus function with period 3.3, with
    Gaussian noise added, and has mean=0 and SD=1."""

    # Since I cannot reproduce the random seed for SciLab.
    # I just grabbed the values from the pdf.
    if default:
        X = [1.0135518, -0.7113242, -0.3906069, 1.565203, 0.0439317,
             -1.1656093, 1.0701692, 1.0825379, -1.2239744, -0.0321446,
             1.1815997, -1.4969448, -0.7455299, 1.0973884, -0.2188716,
             -1.0719573, 0.9922009, 0.4374216, -1.6880219, 0.2609807]
    else:
        N = 20
        t = np.arange(1, N + 1)
        # Sine function with period length 3.3
        X = np.sin(2.0 * np.pi * t / 3.3)
        # Original take a float, but numpy take only integers.
        np.random.seed(123456789)  # 0.123456789

        # Gaussian noise.
        noise = 0.2 * np.random.randn(N)
        X = X + noise
        # Remove mean value.
        X = X - X.mean()
        # Normalize to std=1.
        X = X / X.std()

    return X


# Figure.
X = gen_series()
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 6))
ax.plot(X, '-b.')
ax.set_xlabel('Time series X(t) vs T')
plt.show()

# Some constants.
N = len(X)
M = 4  # Windows size.

# Using a toeplitz matrix
covX = lagcorr(X, X, M=4)  # CovX = [1.0, -0.2867228, -0.6546009, 0.6872825]
C = toeplitz(covX)


# Using shifted time series.
def shift(arr, n, order='forward'):
    if order == 'forward':
        shifted = arr[n:] + [0] * n
    elif order == 'reversed':
        shifted = [0] * n + arr[:-n]
    else:
        print("Order %s not recognized.  Try forward or reversed" % order)

    return shifted

# Embedded Time Series.
Y = np.c_[X, shift(X, 1), shift(X, 2), shift(X, 3)]

# If using the second method.
if 0:
    C = Y.copy()

# Computing the eigenvalues (lambda) eigenvectors (rho) or C:
lamb, rho = np.linalg.eig(C)

"""The columns of the matrix rho are the eigenvectors.  The leading eigenvalue
is in the first row (lambda1 = 2.26); the corresponding leading eigenvector is
in the first column of RHO.  The second eigenvalue is in the second row
(lambda2 = 1.41) and its eigenvector in the second column, etc."""

fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(range(M), rho, marker='.')
ax.legend(['1', '2', '3', '4'], numpoints=1)
ax.set_title(r'Eigenvectors "$\rho$"')
ax.plot(range(M), [0] * len(rho), 'k--')


"""Principal components (The eigenvectors--EOFs) of the matrix C tell us
something about the temporal covariance of the time series, measured at
different lags, as explained above.  We can use them to construct the
principal components of the time series.  The principal components are again
time series, of the same length as the "embedded" time series (i.e., the
matrix Y).  The difference with the embedded time series Y is that we will
choose a different coordinate system to plot each point."""

PC = np.dot(Y, rho)

"""Note that the four columns of the matrix PC are the principal components
PC1, PC2, PC3 and PC4.  They are ordered in the same way as the eigenvectors
are ordered in the matrix W; so the 1st column is PC1, the 2nd column is PC2,
etc."""

"""An important difference between PC and Y is that the columns of PC do not
correspond to different time lags.  Rather, the original values of Y have been
transformed (i.e., projected in a new coordinate system) in order to gather
most of the variance in the first PC, most of the remaining variance in the
second PC, and so on.  Importantly, the PCs are "orthogonal" at lag zero, i.e.,
there is no covariance between the PCs (although there is covariance between
the PCs at nonzero lags).  This is most clearly seen by computing the
variance-covariance matrix of the PC matrix:"""

#np.dot(PC.T, PC / N)

# This is not necessary, but shift takes a list instead of array and I am lazy.
PC0, PC1, PC2, PC3 = map(list, (PC[:, 0], PC[:, 1], PC[:, 2], PC[:, 3]))

# Plotting the PCs.
fig, (ax0, ax1, ax2, ax3) = plt.subplots(nrows=4, ncols=1, figsize=(8, 6))
ax0.set_title(r'Principal components PC vs T')
ax0.plot(range(N), PC0, '.-')
ax1.plot(range(N), PC1, '.-')
ax2.plot(range(N), PC2, '.-')
ax3.plot(range(N), PC3, '.-')

"""Two things are important to notice: (1) the third and fourth PCs contain
very little variance (we already knew that from the eigenvalues); and (2) the
first two PCs are in "quadrature", that is they oscillate as a sine and a
cosine: with the same period but with a 1/4 period phase difference."""

# Reconstruction of the time series.
Z0 = np.c_[PC0, shift(PC0, 1, 'reversed'),
           shift(PC0, 2, 'reversed'),
           shift(PC0, 3, 'reversed')]


Z1 = np.c_[PC1, shift(PC1, 1, 'reversed'),
           shift(PC1, 2, 'reversed'),
           shift(PC1, 3, 'reversed')]

Z2 = np.c_[PC2, shift(PC2, 1, 'reversed'),
           shift(PC2, 2, 'reversed'),
           shift(PC2, 3, 'reversed')]

Z3 = np.c_[PC3, shift(PC3, 1, 'reversed'),
           shift(PC3, 2, 'reversed'),
           shift(PC3, 3, 'reversed')]

RC0 = np.dot(Z0, rho[:, 0]) / M
RC1 = np.dot(Z1, rho[:, 1]) / M
RC2 = np.dot(Z2, rho[:, 2]) / M
RC3 = np.dot(Z3, rho[:, 3]) / M

# Plotting the RCs.
fig, (ax0, ax1, ax2, ax3) = plt.subplots(nrows=4, ncols=1, figsize=(8, 6))
ax0.set_title(r'Reconstructed components RC vs T')
ax0.plot(range(N), RC0, 'r.-')
ax1.plot(range(N), RC1, 'r.-')
ax2.plot(range(N), RC2, 'r.-')
ax3.plot(range(N), RC3, 'r.-')

"""The figure shows that (1) the first two RCs are not in quadrature but in
phase; (2) the first two RCs contain practically all variance of the time
series (we already knew that based on the eigenvalues); the RC4 seems to
describe a "trend" in the data (that might indicate that the used random number
generator is not that good...)."""

# Plotting comparisons.
fig, (ax0, ax1) = plt.subplots(nrows=2, ncols=1, figsize=(8, 6))
ax0.set_title(r'Reconstructed components RC vs T')
ax0.plot(range(N), X, 'b.-', label='Original')
ax0.plot(range(N), RC0 + RC1 + RC2 + RC3, 'r.-', label='RC0-3')
ax0.legend(numpoints=1)
ax1.plot(range(N), X, 'b.-', label='Original')
ax1.plot(range(N), RC0 + RC1, 'r.-', label='RC0-1')
ax1.legend(numpoints=1)

"""The conclusion from these reconstructions is that we have reduced the time
series to oscillatory components (corresponding to the first two eigenvalues)
and two noise components (corresponding to the third and fourth eigenvalues).
In fact, in the second figure, we have used the found RCs to "filter" the time
series: by using less that the total number of RCs, we filter out a part of the
time series that we suppose to be noise rather than signal.  Here, we used the
first two RCs, which we knew define an oscillatory signal due to the phase
quadrature of the corresponding PCs, and their similar eigenvalues (although
this remains to be verified by statistical tests!)."""
