#!/usr/bin/env python
# -*- coding: utf-8 -*-

# import logging, sys
# logging.disable(sys.maxsize)
# import warnings
# warnings.filterwarnings("ignore", category=RuntimeWarning) 

from tqdm import tqdm, trange

import os
import numpy as np
import matplotlib.pyplot as plt

## astylo
from astylo.bio import read_fits, write_fits
from astylo.alib import fixwcs
from astylo.plot import plot2d

## Local
from param import (
	src, Nmc, 
	phot, phot0, path_phot, path_cal, fits_phot, fits_phot0, 
	path_out, path_tmp, path_slices, path_tests, verbose, 
)

def correview(Rx, Ry, yerr, xerr, source=None, threshold=None, \
	title="Correlation", xlabel="x", ylabel="y", figsize=None, \
	savename='correlation.png'):
	
	fig, ax = plt.subplots(figsize=figsize)
	
	x = np.arange(0,10,.01)

	for i, s in enumerate(source):

		## linear fit
		fx = Rx[i].ravel()
		fy = Ry[i].ravel()
		valid = ~(np.isnan(fx) | np.isnan(fy))
		popt, pcov = optimize.curve_fit(f_lin, fx[valid], fy[valid])
		# print("popt = ", popt)
		perr = np.sqrt(np.diag(pcov))
		# print("perr = ", perr)

		## draw nstd-sigma intervals
		nstd = threshold
		popt1 = popt + nstd * perr
		popt2 = popt - nstd * perr
		f1 = f_lin(x, *popt1)
		f2 = f_lin(x, *popt2)
		# ax.fill_between(x, f1, f2, facecolor='grey', alpha=.25, \
			# label="{}-sigma intervals".format(threshold))
		ax.errorbar(x, f_lin(x, *popt), \
			c=COList[i], lw=0.5, ecolor='grey', \
			label="{}: Y = {:.2f} * X + {:.2f}".format(s, *popt))
		if s=='M82':
			ax.errorbar(fx, fy, None, None, \
			fmt='.', markeredgewidth=.01, capsize=2., \
			c=COList[i], ecolor=COList[i], elinewidth=.5, ls=None, \
			label=s)
		## scatter
		else:
			ax.errorbar(fx, fy, yerr[i].ravel(), xerr[i].ravel(), \
			fmt='.', markeredgewidth=.01, capsize=2., \
			c=COList[i], ecolor=COList[i], elinewidth=.5, ls=None, \
			label=s)

	ax.set_title(title)
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	ax.set_xscale('log', nonposx='clip')
	ax.set_yscale('log', nonposy='clip')
	ax.set_xlim(.1, 3.)
	ax.set_ylim(.6, 10.)
	ax.legend(loc='upper left')

	fig.savefig(savename)

def calcuratio(bd1, bd2, bdPATH):

	rationame = bd1+'_'+bd2
	hdr = read_fits(bdPATH+bd1+'/'+bd1+'_0', False)[1]
	Nx = hdr['NAXIS1']
	Ny = hdr['NAXIS2']
	hyperatio = []
	for i in range(N_mc):
		data1 = read_fits(bdPATH+bd1+'/'+bd1+'_'+str(i), False)[0]
		data2 = read_fits(bdPATH+bd2+'/'+bd2+'_'+str(i), False)[0]
		ratio = np.full((Ny, Nx), np.nan)
		for y in range(Ny):
			for x in range(Nx):
				if data2[y,x]!=0 and data1[y,x]!=0:
					ratio[y,x] = data1[y,x] / data2[y,x]
		if i==0:
			ratio0 = ratio
			write_fits(bdPATH+rationame, ratio0, hdr)
		hyperatio.append(ratio)
	hyperatio = np.array(hyperatio)
	error = calunc(hyperatio, [Ny, Nx], bdPATH+rationame+'_unc', hdr)
	
	return ratio0, error


N_mc = 20
threshold = 2.
source = ['M82']
bands = ['I62', 'I77', 'I86', 'I112', 'I127', 'I170']

r1, r2, r3, r4, r5, r6, r7 = [], [], [], [], [], [], []
r1err, r2err, r3err, r4err, r5err ,r6err, r7err = [], [], [], [], [], [], []
for s in source:

	band_path = 'fitting/'+s+'/Intensities/'
	
	## calculate band ratios
	r1label = "$I_{6.2}\,/\,I_{11.3}$"
	r = calcuratio(bands[0], bands[3], band_path)
	r1.append(r[0])
	r1err.append(r[1])

	r2label = "$I_{7.7}\,/\,I_{11.3}$"
	r = calcuratio(bands[1], bands[3], band_path)
	r2.append(r[0])
	r2err.append(r[1])

	r3label = "$I_{8.6}\,/\,I_{11.3}$"
	r = calcuratio(bands[2], bands[3], band_path)
	r3.append(r[0])
	r3err.append(r[1])

	r4label = "$I_{12.7}\,/\,I_{11.3}$"
	r = calcuratio(bands[4], bands[3], band_path)
	r4.append(r[0])
	r4err.append(r[1])

	r5label = "$I_{17.0}\,/\,I_{11.3}$"
	r = calcuratio(bands[5], bands[3], band_path)
	r5.append(r[0])
	r5err.append(r[1])

	r6label = "$I_{7.7}\,/\,I_{6.2}$"
	r = calcuratio(bands[1], bands[0], band_path)
	r6.append(r[0])
	r6err.append(r[1])

	r7label = "$I_{8.6}\,/\,I_{6.2}$"
	r = calcuratio(bands[2], bands[0], band_path)
	r7.append(r[0])
	r7err.append(r[1])

cor_path = 'correlations/'
correview(r1, r2, r2err, r1err, source, 3., \
	None, r1label, r2label, savename=cor_path+'correlation_1.png')
correview(r3, r2, r2err, r3err, source, 3., \
	None, r3label, r2label, savename=cor_path+'correlation_2.png')
correview(r2, r6, r6err, r2err, source, 3., \
	None, r2label, r7label, savename=cor_path+'correlation_3.png')
correview(r2, r7, r7err, r2err, source, 3., \
	None, r2label, r7label, savename=cor_path+'correlation_4.png')
correview(r4, r1, r1err, r4err, source, 3., \
	None, r4label, r1label, savename=cor_path+'correlation_5.png')
correview(r5, r1, r1err, r5err, source, 3., \
	None, r5label, r1label, savename=cor_path+'correlation_6.png')
correview(r5, r6, r6err, r5err, source, 3., \
	None, r5label, r6label, savename=cor_path+'correlation_7.png')

plt.show()
