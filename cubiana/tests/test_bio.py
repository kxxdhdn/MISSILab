#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os
testdir = os.path.dirname(os.path.abspath(__file__))

from astropy.io import fits

## Local
sys.path.insert(0, testdir+'/..') ## astylo path
from astylo.bio import read_fits

## Set path
datdir = testdir+'/dat/'
outdir = testdir+'/out/'


## TEST read_fits
##----------------
ds = read_fits(datdir+'M82_09_SL2')
print(ds.data.shape)

