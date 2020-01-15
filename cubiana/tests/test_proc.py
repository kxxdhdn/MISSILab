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


## TEST stat. dist. (impro.rand_splitnorm)
##-----------------------------------------
# a1 = 1
# a2 = 3
# lam = a1
# tau = a2 / a1

# b = []
# i=1
# for N in range(10000):
# 	theta = np.random.normal(0., 1.) # N(0, 1)
# 	flag = np.random.random() # Uniform dist
# 	if flag<1/(1+tau):
# 		b.append(-abs(theta)*a1)
# 		i+=1
# 	else:
# 		b.append(abs(theta)*a2)
# b = np.array(b)

# print(i/10000)
# print(1/(1+tau))

# hist, bins = np.histogram(b, bins=500)
# center = (bins[:-1] + bins[1:]) / 2
# width = 0.7 * (bins[1] - bins[0])
# plt.bar(center, hist, align='center', width=width)
# plt.show()
# exit()

## TEST sextract + imontage
##--------------------------
path_data = '/Users/dhu/Data/AKARI/data/'
build_tmp = 'out/cubuild/'
if not os.path.exists(build_tmp):
	os.makedirs(build_tmp)

obs_id = ['1420415.1', '1420415.2']
N3 = ['F011213824_N002', 'F011213865_N002']
slit = ['Ns', 'Nh']

par_obs = []
out_tmp = []
for i, obs in enumerate(obs_id):
	for s in slit:
		par_obs.append([obs, s, N3[i]])
		out_tmp.append(build_tmp + '/M83_' + obs + '_' + s)

## Simple/functional tests


unc_out = []
for j, par in enumerate(par_obs):
	sext = sextract(path_data, par)
	cube = sext.spec_build(out_tmp[j], sig_pt=1./3600) # Add pointing accuracy
	# unc_out.append(out_tmp[j]+'_unc') # Symmetric unc
	unc_out.append([out_tmp[j]+'_unc_N', out_tmp[j]+'_unc_P']) # Asymmetric unc
	# print(cube.shape, '\n', sext.wave())

mont = imontage(out_tmp, out_tmp[0])
# print(mont.footprint())
# mont.clean()
# mont.combine('out/M83_IRC', ulist=unc_out, dist='norm') # Symmetric unc
mont.combine('out/M83_IRC', ulist=unc_out, dist='splitnorm') # Asymmetric unc
'''

## Monte-Carlo test
Nmc = 6
for i in range(Nmc+1):
	if i==0:
		wunc = True
		sig_pt = 0.
	else:
		wunc = False
		sig_pt = 1./3600 # Add pointing accuracy

	spec = []
	unc_build = []
	for j, par in enumerate(par_obs):
		sext = sextract(path_data, par)
		cube = sext.spec_build(out_tmp[j]+'_'+str(i), write_unc=wunc, sig_pt=sig_pt)
		spec.append(cube)
		## Asymmetric unc cubes (header without pointing shift)
		if i==0:
			unc_build.append([out_tmp[j]+'_unc_N', out_tmp[j]+'_unc_P'])


## Calculate unc_rep
unc_rep = np.nanstd()

for i in range(Nmc+1):
	mont = imontage(out_tmp[], out_tmp[0])
	mont.combine('out/M83_IRC', method='weighted_avg', \
		ulist=unc_build, dist='splitnorm') # Asymmetric unc
'''

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

t_total = time.time()
print("****** total_time = {:.0f} seconds ******".format(t_total - t0))
