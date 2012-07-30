#!/bin/bash
#
# compile.sh
#
# purpose:
# author:   Filipe P. A. Fernandes
# e-mail:   ocefpaf@gmail
# web:      http://ocefpaf.tiddlyspot.com/
# created:  27-Jul-2012
# modified: Fri 27 Jul 2012 06:48:02 PM BRT
#
# obs:
#

f2py mtm2.f -m mtm2 -h mtm2.pyf
f2py -c mtm2.pyf mtm2.f