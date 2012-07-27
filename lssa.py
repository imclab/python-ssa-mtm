# -*- coding: utf-8 -*-
#
#
# lssa.py
#
# purpose:  Tutorial on lssa
# author:   Filipe P. A. Fernandes
# e-mail:   ocefpaf@gmail
# web:      http://ocefpaf.tiddlyspot.com/
# created:  16-Jul-2012
# modified: Fri 27 Jul 2012 05:32:29 PM BRT
#
# obs:  Least-squares spectral analysis
# http://en.wikipedia.org/wiki/Least-squares_spectral_analysis
#

"""Lomb/Scargle Periodogram
It does not require even spaced data.  Therefore, it is used for sediment and
ice core data."""

r"""Lomb (1976) and Scargle (1982) improve on the simple periodogram by a
slight alteration.  What they showed is that if the cosine and sine
coefficients are normalized separately then the classic periodogram can be used
with unevenly spaced data, and yet the statistical behavior of the power is
identical to the behavior you would expect if you had evenly-spaced points.

To calculate the Lomb-Scargle periodogram of a data set $(t_k, y_k)$ we define,
for every frequency `f`, the time constant $\tau$ by:

$$ \tan(4\pi\tau) = \frac{\sum\sin (4\pi f t_k)}{\sum\cos(4\pi f t_k)} $$,

Then the Lomb-Scargle periodogram estimate of the spectral power $P(f)$ at
frequency $f$ is given by:

$$ P(f) = \frac{1}{2\sigma^2}\left\{\frac{\left[ \sum_k(y_k -
\bar{y})\cos 2\pi f(t_k-\tau) \right]^2}{\sum_k\cos^2 2\pi f(t_k-\tau)} +
\frac{\left[ \sum_k(y_k - \bar{y})\sin 2\pi f(t_k-\tau) \right]^2}
{\sum_k\sin^2 2\pi f(t_k-\tau)}\right\}$$,

This equation is less imposing than it looks.  It has two terms, one for the
cosine transform, the other for the sine transform.  Each term is normalized
separately.  The only complication is that each frequency uses a different
offset $\tau$.  Other than these changes, the equation looks just like ab
ordinary digital Fourier transform.

The Lomb-Scargle method has several advantages over the classic periodogram.
One, obviously, is that paleoclimate data are not evenly spaced.  Although this
can be handle by interpolation, the statistical effects of such interpolation
can be complicated.  Secondly, there is a limit to the ordinary periodogram
that comes about from a process called aliasing.  What this means is that two
signals of different frequencies can have identical sets of values if the
samples are taken at exactly even spacing.  With unevenly-spaced data, this
effect can be substantially reduced.  The net result is that the Lomb-Scargle
periodogram can measure frequencies that would be aliased in evenly-spaced
data.
"""

import numpy as np
import matplotlib.pyplot as plt


def lomb(t, y, freq):
    r"""Calculates Lomb periodogram."""
    # Sets constants.
    nfreq = len(freq)
    fmax, fmin = freq[-1], freq[0]
    power = np.zeros(nfreq)
    f4pi = freq * 4 * np.pi
    pi2 = np.pi * 2.
    n = len(y)
    cosarg = np.zeros(n)
    sinarg = np.zeros(n)
    argu = np.zeros(n)
    var = np.cov(y)  # Variance.
    yn = y - y.mean()
    # Do one Lomb loop.
    for fi in range(nfreq):
        sinsum = np.sum(np.sin(f4pi[fi]) * t)
        cossum = np.sum(np.cos(f4pi[fi]) * t)

        tau = np.arctan2(sinsum, cossum)
        argu = pi2 * freq[fi] * (t - tau)

        cosarg = np.cos(argu)
        cfi = np.sum(yn * cosarg)
        cosnorm = np.sum(cosarg ** 2)

        sinarg = np.sin(argu)
        sfi = np.sum(yn * sinarg)
        sinnorm = np.sum(sinarg ** 2)

        power[fi] = (cfi ** 2 / cosnorm + sfi ** 2 / sinnorm) / 2 * var

    return power

# Tutorial
rand = np.random.rand
age = np.arange(0, 601)
ager = age + 0.3 * rand(age.size) - 0.15
ager[0] = age[0]
ager[600] = age[600]
depth = age / 10.  # Creates depth between 0 and 60.
bkg = np.interp(ager, np.arange(0, 601, 10), rand(61))
# Fake Frequencies at 95 and 127 kyr.
f1, f2 = 1. / 95, 1. / 125
sig = np.cos(2 * np.pi * f1 * ager) + np.cos(2 * np.pi * f2 * ager + np.pi)
o18 = sig + bkg

# Pick frequencies to evaluate spectral power.
freq = np.arange(0, 0.02 + 0.0001, 0.0001)

#from scipy.signal.spectral import lombscargle
#power = lombscargle(age, o18, freq)
power = lomb(age, o18, freq)
power[0] = 0

# Normalize to average 1.
power = power / np.std(power)

# Plot the results.
fig, ax = plt.subplots(nrows=1, ncols=1)
ax.plot(freq, power)
ax.set_title('Lomb tutorial')
ax.set_xlabel(r'Frequencies [cycles kyr$^{-1}$]')
ax.set_ylabel('Spectral power')
