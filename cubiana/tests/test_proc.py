#!/usr/bin/env python
# -*- coding: utf-8 -*-

import time
t0 = time.time()
from tqdm import tqdm, trange

import sys, os, logging
testdir = os.path.dirname(os.path.abspath(__file__))+'/'
# logging.disable(sys.maxsize)

import numpy as np
import matplotlib.pyplot as plt

## Local
sys.path.insert(0, testdir+'/..') ## astylo path
from astylo.bio import read_fits, write_fits
from astylo.astrolib import fixwcs
from astylo.proc import (
	wclean, interfill, hswarp, iswarp, 
	islice, imontage, concatenate, sextract, 
)
from astylo.plot import plot2d

## Set path
datdir = testdir+'dat/'
outdir = testdir+'out/'


## TEST iswarp
##-----------------
## Images
allist = [datdir+'M82_04_SL2', datdir+'M82_08_SL2', datdir+'M82_09_SL2', datdir+'M82_09_SL1']
SL2s = [datdir+'M82_04_SL2', datdir+'M82_08_SL2', datdir+'M82_09_SL2']
## Ref header
# refheader = fixwcs(datdir+'M82_LL2_fp').header

## MC usage
tmpdir = outdir+'tmp/'
swp = iswarp(allist, '9:55:52,69:40:45', '1.67', 
	verbose=False, tmpdir=tmpdir)
Nmc = 2
for j in trange(Nmc+1, leave=True, 
	desc='<iswarp> MC test'):
	if j==0:
		comb = swp.combine(SL2s, 'wgt_avg', keepedge=True, \
			filOUT=outdir+'MC_no', tmpdir=tmpdir+'MC_no/')
	else:
		comb = swp.combine(SL2s, 'wgt_avg', keepedge=True, uncpdf='norm', \
			filOUT=outdir+'MC_'+str(j), tmpdir=tmpdir+'MC_'+str(j)+'/')
	tqdm.write(str(comb.image.shape))

## TEST hswarp
##-------------
# oldimage = read_fits(datdir+'M82_09_SL1').data[0]
# oldheader = fixwcs(datdir+'M82_09_SL1').header
# refheader = fixwcs(datdir+'M82_LL2_fp').header

# sw = hswarp(oldimage, oldheader, refheader, \
# 	keepedge=False, verbose=False, tmpdir=outdir+'tmp_sw/')
# print(sw.image.shape)


## TEST sextract + imontage
##--------------------------
path_data = '/Users/dhu/Data/AKARI/data/'
path_build = outdir+'cubuild/'
if not os.path.exists(path_build):
	os.makedirs(path_build)

obs_id = ['1420415.1', '1420415.2']
N3 = ['F011213824_N002', 'F011213865_N002']
slit = ['Ns', 'Nh']

par_obs = []
out_tmp = []
for i, obs in enumerate(obs_id):
	for s in slit:
		par_obs.append([obs, s, N3[i]])
		out_tmp.append(path_build + 'M83_' + obs + '_' + s)

## Monte-Carlo test
print('sextract: Building cube from slit extraction...')
Nmc = 6

unc_build = []
for j in range(Nmc+1):
	for i, par in enumerate(par_obs):
		sext = sextract(path_data, par)
		if j==0:
			cube = sext.spec_build(out_tmp[i], \
				write_unc=True, sig_pt=0.)
			## Symmetric unc cubes (header without pointing shift)
			# unc_build.append(out_tmp[i]+'_unc')
			## Asymmetric unc cubes (header without pointing shift)
			unc_build.append([out_tmp[i]+'_unc_N', out_tmp[i]+'_unc_P'])
		else:
			## Add pointing accuracy
			cube = sext.spec_build(out_tmp[i]+'_'+str(j), \
				write_unc=True, sig_pt=1./3600)
print('sextract: Building cube from slit extraction...[done]')

ref_irs = outdir+'M83_IRS'
mont = imontage(out_tmp, ref_irs, None, 'ref', 3)
mont.make()

# print(mont.footprint())
# mont.clean()
# mont.combine(outdir+'M83_IRC', method='average', write_mc=True, \
# 	do_rep=True, Nmc=Nmc, filUNC=unc_build, dist='norm') # Symmetric unc
mont.combine(outdir+'M83_IRC', method='wgt_avg', write_mc=True, \
	do_rep=True, Nmc=Nmc, filUNC=unc_build, dist='splitnorm') # Asymmetric unc

## Skip rep
# im0 = mont.combine(outdir+'M83_IRC', method='average', \
# 	filUNC='not_None', do_rep=False)
# im0 = mont.combine(outdir+'M83_IRC', method='wgt_avg', \
# 	filUNC='not_None', do_rep=False)

## Reprendre MC
hyperim = []
for j in range(Nmc):
	out_tmp_mc = []
	build_unc_mc = []
	for f in out_tmp:
		f_mc = f+'_'+str(j+1)
		out_tmp_mc.append(f_mc)
		build_unc_mc.append(f_mc+'_unc')
		# build_unc_mc.append([f_mc+'_unc_N', f_mc+'_unc_P'])
	mont_mc = imontage(out_tmp_mc, None, mont.hdr_ref, 'ref')
	mont_mc.make()
	# hyperim.append(mont_mc.combine(path_build+'M83_IRC_'+str(j+1), write_mc=True, \
	# 		method='average', Nmc=Nmc, filUNC=build_unc_mc, dist='norm'))
	hyperim.append(mont_mc.combine(path_build+'M83_IRC_'+str(j+1), write_mc=True, \
			method='wgt_avg', Nmc=Nmc, filUNC=build_unc_mc, dist='splitnorm'))
	print(read_fits(path_build+'M83_IRC_'+str(j+1)).data.shape)
	hyperim.append(read_fits(path_build+'M83_IRC_'+str(j+1)).data)
hyperim = np.array(hyperim)
unc = np.nanstd(hyperim, axis=0)
write_fits(outdir+'M83_IRC_unc', mont.hdr_ref, unc, mont.wvl)


## TEST concatenate
##------------------

file = [outdir+'M83_IRS']
file.append(outdir+'M83_IRC')
concatenate(file, outdir+'M83')

func = [outdir+'M83_IRS_unc']
func.append(outdir+'M83_IRC_unc')
concatenate(func, outdir+'M83_unc')

ds = read_fits(outdir+'M83')
uds = read_fits(outdir+'M83_unc')
data = ds.data
wvl = ds.wave
unc = uds.data

ma = np.ma.MaskedArray(data, mask=np.isnan(data))
mask_any = ma.mask.any(axis=0)
Ny = mask_any.shape[0]
Nx = mask_any.shape[1]
for y in range(Ny):
	for x in range(Nx):
		# fig, ax = plt.subplots()
		# ax.errorbar(wvl, data[:,y,x])
		# ax.set_xscale('symlog')
		# ax.set_yscale('symlog')
		if mask_any[y,x]==False:
			plot2d(wvl, data[:,y,x], yerr=unc[:,y,x], xlog=1, ylog=1)
plt.show()


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

## TEST FITS ref point shift (impro.crop)
##----------------------------------------
# filIN = 'data/M83_0'
# filOUT = 'out/M83_ref_shift'

# w = fixwcs(filIN).WCS
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
