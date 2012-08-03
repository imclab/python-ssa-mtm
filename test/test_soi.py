# -*- coding: utf-8 -*-
#
#
# test_soi.py
#
# purpose:  Test SSA-MTM
# author:   Filipe P. A. Fernandes
# e-mail:   ocefpaf@gmail
# web:      http://ocefpaf.tiddlyspot.com/
# created:  03-Aug-2012
# modified: Fri 03 Aug 2012 06:54:39 PM BRT
#
# obs: Download latest SOI data and try to reproduce Ghil et al. 2002 figures.
#

import urllib2
import numpy as np
from datetime import date
from StringIO import StringIO


def download_soi(url="http://www.cpc.ncep.noaa.gov/data/indices/soi"):
    r"""Download and clean SOI data from NOAA."""
    request = urllib2.Request(url)
    data = urllib2.urlopen(request).read()
    s = StringIO(data)

    # Current year.
    start, stop = [], []
    year_start = '1951'
    year_end = str(date.today().timetuple()[0])
    for k, line in enumerate(s.readlines()):
        if year_start in line:
            start.append(k)
        elif year_end in line:
            stop.append(k)
        else:
            pass

    delimiter = [4, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6]
    dtype = "i8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8"

    # First pass.
    s.seek(0)
    anomaly = np.genfromtxt(s, delimiter=delimiter,
                            dtype=dtype,
                            autostrip=True,
                            skip_header=start[0] - 1,  # Include the header.
                            skip_footer=k - stop[0],
                            #skip_footer=stop[1] - stop[0],
                            missing_values='-999.9',  # FIXME
                            filling_values=np.NaN,
                            names=True)

    # Second pass.
    s.seek(0)
    standardized = np.genfromtxt(s, delimiter=delimiter,
                                 dtype=dtype,
                                 autostrip=True,
                                 skip_header=start[1] - 1,  # Inc. the header.
                                 skip_footer=k - stop[1],
                                 missing_values='-999.9',  # FIXME
                                 filling_values=np.NaN,
                                 names=True)
    return anomaly, standardized
