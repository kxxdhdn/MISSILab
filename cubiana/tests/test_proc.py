#!/usr/bin/env python
# -*- coding: utf-8 -*-

import time
t0 = time.time()

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+'/..')

import numpy as np
import matplotlib.pyplot as plt

from astylo.bio import read_fits, write_fits, ext_wcs
from astylo.proc import wclean, interfill, sextract, imontage


## TEST sextract + imontage
##--------------------------
path_data = '/Users/dhu/Data/AKARI/data/'

obs_id = ['1420415.1', '1420415.2']
N3 = ['F011213824_N002', 'F011213865_N002']
slit = ['Ns', 'Nh']

par_obs = []
file_out = []
for i, obs in enumerate(obs_id):
	for s in slit:
		par_obs.append([obs, N3[i], s])
		file_out.append('out/M83_' + obs + '_' + s)

spec = []
unc_out = []
for j, par in enumerate(par_obs):
	spec.append(sextract(file_out[j], path_data, par))
	unc_out.append(file_out[j]+'_unc')
# print(spec[0].image().shape, spec[0].wave())

# mont = imontage(file_out, file_out[0])
mont = imontage(file_out[0:2], file_out[0])
mont.combine(filOUT='out/M83_IRC', ulist=unc_out)

## TEST FITS ref point shift (impro.crop)
##----------------------------------------
# filIN = 'data/M83_0'
# filOUT = 'out/M83_ref_shift'

# w = ext_wcs(filIN).WCS
# hdr, data, wave = read_fits(filIN)
# hdr['CRPIX1'] = hdr['NAXIS1'] / 2.
# hdr['CRPIX2'] = hdr['NAXIS2'] / 2.
# pix = np.array([(hdr['CRPIX1'], hdr['CRPIX2'])])
# hdr['CRVAL1'], hdr['CRVAL2'] = w.all_pix2world(pix, 1)[0]

# write_fits(filOUT, hdr, data, wave)

## TEST interfill
##----------------
# Nmc = 2

# filIN = 'data/IC10_SL1'
# filOUT = 'out/IC10_SL1'

# hdr, im, wvl = read_fits(filIN, wmod=1)
# unc = read_fits(filIN+'_unc', wmod=1)[1]

# mu, sigma = 0., 1.
# hyperim = []
# for i in range(Nmc+1):
# 	if i==0:
# 		newim = interfill(im, axis=1)
# 		write_fits(filOUT, hdr, newim, wvl)
# 	else:
# 		iunc = im + unc * np.random.normal(mu, sigma, im.shape)
# 		hyperim.append(interfill(iunc, axis=1))
# hyperim = np.array(hyperim)
# print(hyperim.shape)
# newunc = MCerror(hyperim)
# write_fits(filOUT+'_unc', hdr, newunc, wvl)


## TEST wclean
##-------------
# path = '/Users/dhu/data/pahpedia/M83/output/'
# file = 'M83_0'

# wclean(path+file, cmod='closest_right', filOUT=path+file+'_wclean_test')

# t_total = time.time()
# print("****** total_time = {:.0f} seconds ******".format(t_total - t0))
