#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

SOURCE: VV114 (IC 1623, ARP 236)

"""

import os

## rapyuta
from rapyuta.astrom import fixwcs
from rapyuta.inout import fitsext
from rapyuta.imaging import iswarp

## Local
from buildinfo import ( src, Nmc, verbose,
                        chnl, sub_SL, sub_LL, sub_SH, sub_LH,
                        path_irs, path_tmp, path_out
)

## IRC atlas
coadd_footprint = fixwcs(path_out+src+'_footprint'+fitsext).header

## IRS
# refheader = fixwcs(path_irs+src+'_IRC'+fitsext).header
# swp = iswarp(sum(fits_irs, []), refheader=refheader,
#              tmpdir=path_tmp, verbose=verbose)
# coadd_footprint = swp.refheader
