#!/usr/bin/env python
# -*- coding: utf-8 -*-

import time
t0 = time.time()

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+'/..')

import numpy as np

from astylo.lib import nanavg, nanstd

## TEST nanavg & nanstd
##----------------------
a = np.arange(24).reshape((4,3,2))
data = a*1.
unc = a*.1
wgt = unc/np.sum(unc, axis=0) # Nomalized
# print(data)
# print(unc)
# print(wgt)
# print('--------')

# av1 = np.average(data, axis=0)
# av2 = np.average(data, axis=0, weights=unc, returned=True)
# av3 = np.average(data, axis=0, weights=wgt, returned=True)
# print(av1)
# print(av2)
# print(av3)
# print('--------')

data[0,0,0] = np.nan
unc[0,0,0] = np.nan
data[:,1,1] = np.nan
unc[:,1,1] = np.nan
unc[0,2,1] = 0
# print(data[:,0,0])
# print(data[:,1,1])
# print(data[:,2,1])
# print('--------')

# av1 = nanavg(data, axis=0)
# av2 = nanavg(data, axis=0, weights=unc)
# wgt = unc/np.nansum(unc, axis=0) # Nomalized
# av3 = nanavg(data, axis=0, weights=unc)
# print(av1)
# print(av2)
# print(av3)
# print('--------')

inv_var = 1./unc**2 / np.nansum(1./unc**2, axis=0)
# print(inv_var) # Nomalized
# print(1./unc**2) # Not nomalized
# print(nanavg(data, axis=0, weights=inv_var))
# print(nanavg(data, axis=0, weights=1./unc**2)) # Smallest weighted mean
# print(np.nanmean(data, axis=0))

print(nanstd(data, axis=0, weights=unc))
print(nanstd(data, axis=0, weights=1./unc**2))
print(np.nanstd(data, axis=0))


t_total = time.time()
print("****** total_time = {:.0f} seconds ******".format(t_total - t0))
