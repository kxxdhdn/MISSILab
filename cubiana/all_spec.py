#!/usr/bin/env python
# -*- coding: utf-8 -*-

# import logging, sys
# logging.disable(sys.maxsize)
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 

from tqdm import tqdm, trange

import os
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import gmean
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, LogFormatter, NullFormatter
## astylo
from astylo.iolib import fclean, read_fits, write_fits
from astylo.ipro import hextract, hswarp, iconvolve, iswarp, concatenate
from astylo.calib import intercalib
from astylo.mlib import f_lin, f_lin0
from astylo.alib import pix2sr, get_pc, fixwcs

## Local
from param import (
	src, Nmc, path_idl, path_ker, fits_ker, csv_ker, 
	phot, phot0, path_phot, path_cal, fits_phot, fits_phot0, 
	path_out, path_tmp, path_slices, path_tests, verbose, 
)

##---------------------------
##       Initialisation
##---------------------------
Nmc = 2
'''
##---------------------------
##     Rescale pixel size
##---------------------------
ds0 = read_fits(path_out+src+'_IRC_0')
oldimage = ds0.data
wvl = ds0.wave
oldheader = fixwcs(path_out+src+'_IRC_0').header
write_fits(path_tmp+src+'_IRC_fp', oldheader, oldimage[0])

swp = iswarp([path_tmp+src+'_IRC_fp'], \
	pixscale=3., tmpdir=path_tmp, \
	verbose=False)
swp.footprint(path_tmp+'footprint')

for j in trange(Nmc+1, #leave=False, 
	desc='Rescaling pixel size [MC]'):
	## IRC
	swp.combine([path_out+src+'_IRC_'+str(j)], \
		keepedge=True, tmpdir=path_slices, \
		filOUT=path_tmp+src+'_IRC_'+str(j))
	## IRS
	swp.combine([path_out+src+'_IRS_'+str(j)], \
		keepedge=True, tmpdir=path_slices, \
		filOUT=path_tmp+src+'_IRS_'+str(j))

## Cal unc (instr)
##-----------------
mcimage = []
for j in trange(Nmc+1, #leave=False, 
	desc='IRC Reading [MC]'):
	if j==0:
		hd0 = read_fits(path_tmp+src+'_IRC_0')
		header = hd0.header
		wvl = hd0.wave
	else:
		hd = read_fits(path_tmp+src+'_IRC_'+str(j))
		mcimage.append(hd.data)
if Nmc>1:
	mcimage = np.array(mcimage)
	unc = np.nanstd(mcimage, axis=0)
	write_fits(path_tmp+src+'_IRC_unc', header, unc, wvl)

mcimage = []
for j in trange(Nmc+1, #leave=False, 
	desc='IRS Reading [MC]'):
	if j==0:
		hd0 = read_fits(path_tmp+src+'_IRS_0')
		header = hd0.header
		wvl = hd0.wave
	else:
		hd = read_fits(path_tmp+src+'_IRS_'+str(j))
		mcimage.append(hd.data)
if Nmc>1:
	mcimage = np.array(mcimage)
	unc = np.nanstd(mcimage, axis=0)
	write_fits(path_tmp+src+'_IRS_unc', header, unc, wvl)

##---------------------------
##       Concatenate ALL
##---------------------------
for j in trange(Nmc+1, #leave=False, 
	desc='ALL concatenation [MC]'):
	files = [path_tmp+src+'_IRC_'+str(j), path_tmp+src+'_IRS_'+str(j)]
	if j==0:
		func = [path_tmp+src+'_IRC_unc', path_tmp+src+'_IRS_unc']

		concatenate(files, path_tmp+src+'_ma', sort_wave=True)
		data0 = read_fits(path_tmp+src+'_ma').data
		ma = np.ma.array(data0, mask=np.isnan(data0))
		mask_any = ma.mask.any(axis=0)

		concatenate(files, path_out+src+'_0', sort_wave=True)
		ds = read_fits(path_out+src+'_0')
		data = ds.data
		header = ds.header
		wvl = ds.wave
		for k in range(len(wvl)):
			data[k][mask_any] = np.nan
		write_fits(path_out+src+'_0', header, data, wvl)
	else:
		concatenate(files, path_out+src+'_'+str(j), sort_wave=True)
		data = read_fits(path_out+src+'_'+str(j)).data
		for k in range(len(wvl)):
			data[k][mask_any] = np.nan
		write_fits(path_out+src+'_'+str(j), header, data, wvl)

## Cal unc (ALL)
##---------------
mcimage = []
for j in trange(Nmc+1, #leave=False, 
	desc='ALL Reading [MC]'):
	if j==0:
		hd0 = read_fits(path_out+src+'_0')
		header = hd0.header
		wvl = hd0.wave
	else:
		hd = read_fits(path_out+src+'_'+str(j))
		mcimage.append(hd.data)
if Nmc>1:
	mcimage = np.array(mcimage)
	unc = np.nanstd(mcimage, axis=0)
	write_fits(path_out+src+'_unc', header, unc, wvl)
'''

##---------------------------
##       Plot spectra
##---------------------------
# ds = read_fits(path_out+src+'_0', path_out+src+'_unc')
# data = ds.data
# data[:259,:,:] = data[:259,:,:]*data[259,:,:]/data[258,:,:]
# wvl = ds.wave

# Nw = len(wvl)
# unc = ds.unc
# newimage = []
# newwvl = []
# newunc = []
# for k in range(Nw):
# 	if k>6 and k<259:
# 		newimage.append(data[k,:,:])
# 		newunc.append(unc[k,:,:])
# 		newwvl.append(wvl[k])
# 	if k>263 and k<466:
# 		newimage.append(data[k,:,:])
# 		newunc.append(unc[k,:,:])
# 		newwvl.append(wvl[k])
# 	if k>481:
# 		newimage.append(data[k,:,:])
# 		newunc.append(unc[k,:,:])
# 		newwvl.append(wvl[k])
# newimage = np.array(newimage)
# newunc = np.array(newunc)
# newwvl = np.array(newwvl)

# ma = np.ma.array(data, mask=np.isnan(data))
# mask_any = ma.mask.any(axis=0)

# Ny = mask_any.shape[0]
# Nx = mask_any.shape[1]
# for y in range(Ny):
# 	for x in range(Nx):
# 		fig, ax = plt.subplots()
# 		ax.errorbar(wvl, data[:,y,x])
# 		ax.set_xscale('symlog')
# 		ax.set_yscale('symlog')
# 		if mask_any[y,x]==False:
# 			print(data[100,y,x]/unc[100,y,x])
			# plot2d(wvl, data[:,y,x], yerr=unc[:,y,x], xlog=1, ylog=1)

# x = [30,34,38,39,39,41,41,43,43,45,47,49,52,]
# y = [34,29,39,26,32,26,28,30,35,31,26,28,30,]
# image = []
# for i in range(len(x)):
# 	# image.append(newimage[:,y[i],x[i]])
# 	image.append(data[:,y[i],x[i]])

# 	# plot2d(newwvl, newimage[:,y[i],x[i]], yerr=newunc[:,y[i],x[i]], xlog=1, ylog=1)
# # plot2d_m([newwvl]*13, image, xlog=1, ylog=1)
# plot2d_m([wvl]*13, image, xlog=1, ylog=1)
# plt.show()



##---------------------------
##     Inter-Calibration
##---------------------------
## Reproject IRAC
phot = ['IRAC1', 'IRAC2', 'IRAC3', 'IRAC4']
refheader = read_fits(path_out+'footprint').header
F_phot = []
for p in phot:
	ds = read_fits(path_phot+src+'_'+p)
	oldimage = ds.data
	oldheader = ds.header
	swp = hswarp(oldimage, oldheader, refheader, keepedge=False)
	newimage = swp.image
	F_phot.append(newimage)
	write_fits(path_cal+src+'_'+p, refheader, newimage)

## Synthetic photometry
calib = intercalib(path_out+src)
wcen, F_syn, sig = calib.synthetic_photometry(('IRAC1', 'IRAC2', 'IRAC3', 'IRAC4'))
for i in range(4):
	write_fits(path_cal+src+'_sp_'+phot[i], refheader, F_syn[i])

i = 0
syn = read_fits(path_cal+src+'_sp_'+phot[i]).data.reshape((-1,))
pho = read_fits(path_cal+src+'_'+phot[i]).data.reshape((-1,))

mask = ~np.ma.array(syn, 
	mask=np.logical_or(np.isnan(syn), np.isnan(pho))).mask
plt.scatter(syn[mask], pho[mask])
plt.xscale('symlog')
plt.yscale('symlog')
# plt.xlim((0,2.e2))
# plt.ylim((0,2.e2))
plt.show()


