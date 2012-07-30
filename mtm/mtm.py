# -*- coding: utf-8 -*-
#
#
# mtm.py
#
# purpose:  Wrap mtm2.f
# author:   Filipe P. A. Fernandes
# e-mail:   ocefpaf@gmail
# web:      http://ocefpaf.tiddlyspot.com/
# created:  27-Jul-2012
# modified: Fri 27 Jul 2012 06:43:45 PM BRT
#
# obs:  MultiTaper Spectral Analysis.
#

import mtm2
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.font_manager as fm
#from matplotlib import rcParams
#rcParams['font.family'] = 'sans-serif'
#rcParams['font.sans-serif'] = ['Tahoma']

def compute_mtm(arr, dt=1.0, npi=2, nwin=3, f1=0.0, f2=0.0, inorm=0, ispec=1,
                iresh=1, ithresh=3, inoise=0, ilog=1, ismooth=1, isignal=0,
                irecon=0, nsignals=0, iplotresh=1, iplotftest=1, iplotsmoo=1,
                iplotraw=1, iplotconf=1, icon=0):
    r"""This function compute mtm spectrum on a time-series array.  It returns
    five arrays:
    * 3 columns array spec_raw with the frequency, the raw spectrum and ftest
      in the last one.
    * 3 columns array spec_resh with the frequency, the reshape signal
      spectrum and the harmonic.
    * 2 columns array spec_smoo with the frequency and the smooth spectrum.
    * 5 columns array spec_conf with the frequency, the confidence for the
      spectrum obtained.
    * 3 columns array recon with the frequency, the reconstructed signal.

    It has a large  number of parameters which are initially set to some
    values:

    The time step parameter: `dt` = 1.0 in years
        To compute mtm on monthly data, set dt = 0.0833333333 (1/12 year)

    The resolution in `p` multiples of Rayleigh frequency: npi = 2
    The number of tapers: suggested nwin = 2 * npi - 1
        if nwin > (2 * npi - 1) then nwin = 2 * npi - 1

    The spectrum:
        f1 = 0.0  (fmin)
        f2 = fny = 0.5/dt  (fmax)

    The normalization: inorm = 0
        inorm = (0) none (1) N  (2) 1/dt

    The spectrum resolution: ispec = 1
        ispec = (0) high-resolution or (1) adaptive estimate

    Reshaping flag: iresh = 1
        iresh = reshaping (yes=1)?
        if iresh = 1
            itresh = sig. level for harmonic peak detection/reshaping
                (0) 50% (1) 90% (2) 95% (3) 99% (4) 99.5% (5) 99.9%
            ithresh = 3

    Some hypothesis parameters:
        inoise = (0) red noise, (1) white noise, (2) locally-white noise
        inoise = 0
        ilog = red noise background fit: (0) linear,(1) log
        ilog = 1
        ismooth = background estimation: (0) raw, (1) robust'
        ismooth = 1

    Signal assumption:
        isignal = (0) harmonic or narrowband, (1) narrowband, (2) harmonic
        isignal = 0
        if isignal=2 then inoise=2, iplotsmoo=0, iplotconf=0
        if isignal=1 then iplotftest=0 and iplotresh=0

    Reconstruction:
        irecon = 0
        signals to reconstruct: nsignals=0

    Display option parameters:
        (1) raw spectrum                - iplotraw
        (2) reshaped & harmonic spectra - iplotresh
        (3) median smoothed background  - iplotsmoo
        (4) 50% 90% 95% 99% conf levels - iplotconf
        (5) F-test spectrum             - iplotftest
        iplotresh = 1
        iplotftest = 1
        iplotsmoo = 1
        iplotraw = 1
        iplotconf = 1

    Constraint Option:
        icon = (0) min misfit, (1) min norm, (2) min slope, (3) min rough
        icon=0
    """

    if nwin > (2 * npi - 1):
        nwin = 2 * npi - 1
    if f2 == 0:
        f2 = 0.5 / dt
    if isignal == 2:
        inoise = 2
        iplotsmoo = 0
        iplotconf = 0
    if isignal == 1:
        iplotftest = 0
        iplotresh = 0

    spec_raw, spec_resh, spec_smoo, spec_conf, recon_sig = mtm2.mtm_mann(
        arr, dt, npi, nwin, f1, f2, inorm, ispec, iresh, ithresh, inoise, ilog,
        ismooth, isignal, irecon, nsignals, iplotresh, iplotftest, iplotsmoo,
        iplotraw, iplotconf, icon)

    spec_raw = resize_spec(spec_raw)
    spec_resh = resize_spec(spec_resh)
    spec_smoo = resize_spec(spec_smoo)
    spec_conf = resize_spec(spec_conf)
    recon_sig = resize_spec(recon_sig)

    return spec_raw, spec_resh, spec_smoo, spec_conf, recon_sig


def plot_spec(spec_raw):
    fig, ax = plt.subplots()
    #prop = fm.FontProperties(fname='/Library/Fonts/Tahoma.ttf')
    #ax.set_xlabel('Frequency [years]', fontsize=10, color='black',
    #fontproperties=prop)
    m = np.array([1.0, 2.0, 5.0, 10.0, 20.0, 100.0])
    ohm = 1.0 / m
    ax.loglog(spec_raw[:, 0], spec_raw[:, 1], linewidth=0.5)
    for i in range(0, 6):
        ax.text(ohm[i], 200000, str(int(m[i])))
    ax.set_xlabel('Frequency [1/years]', fontsize=12, color='black')
    ax.set_ylabel('Spectral Power', fontsize=12, color='black')
    ax.set_title('MTM Spectrum', fontsize=12, color='black')
    ax.text(2, 200000, 'Period [year]')
    ax.text(.01, 0.01, r'Nino3 from Speedy')
    #ax.axis([0.005, 10, 0.01, 1000000])
    plt.grid(True)
    plt.show()


def resize_spec(arr):
    arr = np.asanyarray(arr)
    if arr.ndim != 2:
        raise ValueError("Array must be 2D.")
    size = 1
    nrow, ncol = arr.shape
    for i in range(1, nrow):
        if arr[i, 0] != 0:
            size = size + 1
    if size != 1:
        arr = arr[0:size, 0:ncol]
    return arr


def find_inequalities(arr, fmin, fmax):
    if3 = 0
    nrow, ncol = arr.shape
    if7 = nrow - 1
    while arr[if3, 0] < fmin and if3 < nrow:
        if3 = if3 + 1
    while arr[if7, 0] > fmax and if7 > 0:
        if7 = if7 - 1
    return if3, if7


if __name__ == '__main__':
    X = [1.0135518, -0.7113242, -0.3906069, 1.565203, 0.0439317, -1.1656093,
         1.0701692, 1.0825379, -1.2239744, -0.0321446, 1.1815997, -1.4969448,
         -0.7455299, 1.0973884, -0.2188716, -1.0719573, 0.9922009, 0.4374216,
         -1.6880219, 0.2609807]

    X = np.atleast_2d(X)
    spec_raw, spec_resh, spec_smoo, spec_conf, recon_sig = compute_mtm(X)
    plot_spec(spec_raw)