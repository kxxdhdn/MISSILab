#!/usr/bin/env python
# -*- coding: utf-8 -*-

import time
t0 = time.time()

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+'/..')

import numpy as np
import matplotlib.pyplot as plt

from astylo.processim import wclean


## TEST wclean
##-------------
path = '/Users/dhu/data/pahpedia/M83/output/'
file = 'M83_0'

wclean(path+file, filOUT=path+file+'_wclean_test')

t_total = time.time()
print("****** total_time = {:.0f} seconds ******".format(t_total - t0))
