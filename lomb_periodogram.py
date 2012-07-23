# -*- coding: utf-8 -*-
#
# spectral_tests.py
#
# purpose:  Lomb spectra analysis test
# author:   Filipe P. A. Fernandes
# e-mail:   ocefpaf@gmail
# web:      http://ocefpaf.tiddlyspot.com/
# created:  04-Jul-2012
# modified: Mon 23 Jul 2012 01:11:04 PM BRT
#
# obs:
#

"""As principais diferenças entre esse tipo de análise e uma análise espectral
tradicional é a natureza de mínimos quadrados para resolver senos e cossenos
separadamente.  Ela é considerada "qualitativa" e não quantitativa.  Mesmo os
picos que não são prováveis, podem e devem ser discutidos na luz de
conhecimento pretérito.  Como um ciclo bem conhecido por exemplo.

Apesar de ser qualitativa é mais recomendado que interpolação quando a série é
heterogênea.  Pois a interpolação vai introduzir muito mais erros e/ou ser
impossível.  Como nesse caso, onde os intervalos varia entre 6 anos e +3mil
anos.

Leia mais sobre aqui:
    http://en.wikipedia.org/wiki/Least-squares_spectral_analysis

Já os dados da Amanda são entre 29 e 30 anos de intervalo, interpolar nesse
caso é seguro.  Porém, a série é curta e os intervalos de confiança largos.
Para se ter uma análise melhor tem-se que tentar a Singular Spectrum Analysis
(SSA).  Já fiz tutoriais para essas, mas não testei com do dados dela."""

import numpy as np
import matplotlib.pyplot as plt


def lombscargle(ages, signal, ofac=4, hifac=1):
    r"""Calculates Lomb-Scargle Periodogram.

    Enter `signal` at times `ages` to compute the periodogram, with
    oversampling factor `ofac` and up to frequencies `hifac` * Nyquist.

    Return frequencies considered `freq`, the associated spectral `power`,
    and estimated significance of the power values `prob`.

    Note: the significance returned is the false alarm probability of the null
    hypothesis, i.e. that the data is composed of independent Gaussian random
    variables.  Low probability values indicate a high degree of significance
    in the associated periodic signal."""

    N, T = len(signal), ages.ptp()

    # Mean and variance.
    mu, s2 = signal.mean(), signal.var()

    # Calculate sampling frequencies.
    start = 1.0 / (T * ofac)
    stop = hifac * N / (2.0 * T)
    dt = 1.0 / (T * ofac)  # Interval for the frequencies.  Can be tweaked.
    freq = np.arange(start, stop + dt, dt)

    # Angular frequencies and constant offsets.
    w = 2.0 * np.pi * freq
    dot = np.dot(w[:, None], ages[None, :])
    A = np.sum(np.sin(2.0 * dot), axis=1)
    B = np.sum(np.cos(2.0 * dot), axis=1)
    tau = np.arctan2(A, B) / (2.0 * w)

    # Spectral power.
    cterm = np.cos(dot - (w * tau)[:, None])
    sterm = np.sin(dot - (w * tau)[:, None])

    ry = (np.sum(np.dot(cterm, np.diag(signal - mu)), axis=1) ** 2.0 /
          np.sum(cterm ** 2, axis=1))
    iy = (np.sum(np.dot(sterm, np.diag(signal - mu)), axis=1) ** 2.0 /
          np.sum(sterm ** 2, axis=1))

    # TODO: Phase (untested!)
    phLS = np.arctan2(ry, iy)

    power = (np.sum(np.dot(cterm, np.diag(signal - mu)), axis=1) ** 2.0 /
             np.sum(cterm ** 2, axis=1) +
             np.sum(np.dot(sterm, np.diag(signal - mu)), axis=1) ** 2.0 /
             np.sum(sterm ** 2, axis=1))

    power /= (2.0 * s2)

    # Estimate of the number of independent frequencies.
    M = 2.0 * len(freq) / ofac

    # Statistical significant of power.
    prob = M * np.exp(-power)
    inds = prob > 0.01
    prob[inds] = 1.0 - (1.0 - np.exp(-power[inds])) ** M

    return freq, power, prob, phLS


# TODO
def revert_signal(ages, phases, ofac=4, hifac=1):
    r"""Reconstruct the signal with a specified frequency."""

    # Same constants.
    max_t, min_t = ages.max(), ages.min()
    dif_t = ages.ptp()
    ave_t = 0.5 * (max_t + min_t)

    # Calculate sampling frequencies.
    N = len(ages)
    start = 1.0 / (dif_t * ofac)
    stop = hifac * N / (2.0 * dif_t)
    dt = 1.0 / (dif_t * ofac)
    freq = np.arange(start, stop + dt, dt)

    # Phase shift with respect to 0.
    arg0 = 2.0 * np.pi * (ave_t + ages[0]) * freq + wtau

    # Phase shift for FFT reconstruction.
    arg1 = 2.0 * np.pi * ave_t * freq + wtau
    ph0 = np.mod(phLS + arg0, 2.0 * np.pi)  # Phase with respect to 0.
    ph1 = np.mod(phLS + arg1, 2.0 * np.pi)  # Phase for complex FFT spectrum.
    return ph0, ph1


def explore_periods():
    print("Click on the desired peak to check its period in Years.")
    x, y = plt.ginput(1)[0]
    period = 1. / x * 1e3
    print("Period clicked is: ~ %s years." % np.int_(period))
