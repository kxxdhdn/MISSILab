#!/usr/bin/env python
# -*- coding: utf-8 -*-

import time
t0 = time.time()

import os
import numpy as np
import matplotlib.pyplot as plt
## astylo
from astylo.bio import read_fits, write_fits, read_ascii
from astylo.proc import islice, icrop, iconvolve, imontage, wclean
from astylo.calib import intercalib, phot2phot
from astylo.lib import fclean, pix2sr
from astylo.plot import plot2d

##---------------------------
##       Initialisation
##---------------------------

Nmc = 6

src = 'M83'
## Not useful if need to modify IDL/conv_prog.pro
# src = input("Input source name: ")

verbose = False

##---------------------------
##         Path Store
##---------------------------

