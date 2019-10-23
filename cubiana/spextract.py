#!/usr/bin/env python
# -*- coding: utf-8 -*-

import time
t0 = time.time()

import os
import numpy as np
## astylo
from astylo.sinout import read_fits, write_fits, WCSextract, read_ascii
from astylo.processim import slicube, crop, iconvolve, project
from astylo.myfunclib import fclean
from astylo.mc import calunc

##---------------------------
##       Initialisation
##---------------------------

Nmc = 2

src = 'M83'
wmod = 1
## Not useful if need to modify IDL/conv_prog.pro
# src = input("Input source name: ")
if src=='IC10':
	wmod = 0

instr = ['IRS', 'IRC']
chnl = ['SL2', 'SL1', 'LL2']
# chnl = ['SL2', 'SL1', 'LL2', 'LL1']

##---------------------------
##         Path Store
##---------------------------

path_data = '/Users/dhu/data/pahpedia/'+src+'/' ###
path_tmp = '/Users/dhu/data/pahpedia/tmp/' ### See also IDL/convolve_image.pro
if not os.path.exists(path_tmp):
	os.makedirs(path_tmp)
path_test = '/Users/dhu/data/pahpedia/tests/' ###
path_ker = '/Users/dhu/data/kernels/' ###
path_out = path_data+'output/' ###vvv See also
if not os.path.exists(path_out):
	os.makedirs(path_out)

## Current dir
path_root = os.getcwd()+'/'
path_idl = path_root+'/IDL/'
## Default files
chnl_ref = path_data+src+'_'+chnl[-1] ###
file_ref = path_data+src+'_ref'
file_all = path_out+src

list_ker = path_out+'kernels_'+src ### See also IDL/conv_prog.pro

##---------------------------
##         Read Data
##---------------------------

## Image
coord = read_ascii(path_root+'data/coord')
for c in coord:
	if c[0]==src:
		ra, dec = float(c[1]), float(c[2])
dx, dy = 15, 15
## Kernel
psf = [2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.]
kernelist = []
for p in psf:
	kernelist.append(path_ker+'Kernel_HiRes_Gauss_0'+ \
		str(p)+'_to_Gauss_06.0') ###

##---------------------------
##         Processing
##---------------------------

## Cropped ref image
##-------------------
cr = crop(filIN=chnl_ref, cen=(ra,dec), size=(dx,dy), \
	wmod=wmod, filOUT=file_ref)
hdr = read_fits(file_ref)[0]
NAXIS1 = hdr['NAXIS1']
NAXIS2 = hdr['NAXIS2']

t_init = time.time()
print("****** init_time = {:.0f} seconds ******".format(t_init - t0))

b0 = input("(Re)do homegeneisation? [y/n] ")
# b0 = 'y'
if b0=='y':
	b1 = input("(Re)do Monte-Carlo? [y/n] ")
	
	## Time offset
	t1 = time.time()

	# b1 = 'y'
	if b1=='y':
		hypercube=[]
		Nmc += 1
	else:
		Nmc = 1
	
	for j in range(Nmc):
		cubi = []
		wavALL = []
		slices = []
		for i, ch in enumerate(chnl):

			file_data = path_data+src+'_'+ch
			file_unc = file_data+'_unc'
			file_slice = path_tmp+src+'_'+ch
			file_out = path_out+src+'_'+ch
			file_conv = file_out+'_conv'

			## Slice cube
			##------------
			# if j==0: # unc not added
			# 	sl = slicube(file_data, file_slice, wmod=wmod)
			# 	wvl = sl.wave()
			# 	slist = sl.slice_names()
			# else: # add unc
			# 	slicube(filIN=file_data, filSL=file_slice, \
			# 		wmod=wmod, uncIN=file_unc, wmod_unc=wmod)

			## Smooth images [IDL]
			##---------------------
			if j==0: # unc not added
				conv = iconvolve(filIN=file_data, filKER=kernelist, \
					saveKER=list_ker, wmod=wmod, \
					filTMP=file_slice, filOUT=file_conv)
				wavALL.extend(conv.wave())
			else: # add unc
				conv = iconvolve(filIN=file_data, filKER=kernelist, \
					saveKER=list_ker, wmod=wmod, uncIN=file_unc, wmod_unc=wmod, \
					filTMP=file_slice, filOUT=file_conv)

			
			conv.do_conv(ipath=path_idl)

			## Reproject convolved cube
			##--------------------------
			pr = project(filIN=file_conv, filREF=file_ref, \
				filTMP=file_slice, filOUT=file_out+'_'+str(j))
			
			slices.extend(pr.slice_names())
			cubi.append(pr.image()) # Here, pr.image is a 3D cube
		
		cube = np.concatenate(cubi)
		NAXIS3 = np.size(wavALL)

		comment = "Homegeneized cube produced by [SPEXTRACT] routine. "

		write_fits(file_all+'_'+str(j), hdr, cube, wave=wavALL, COMMENT=comment)

		if j==0:
			## no uncertainty added cube
			cube0 = cube

			t_no_mc = time.time()
			print("****** no_mc_time = {:.0f} seconds ******".format(t_no_mc - t1))
		else:
			## append hypercube
			hypercube.append(cube)

			## ieme iteration finished
			print("---------------- {} ----------------".format(j))

	if b1=='y':
		## calculate MC uncertainty
		t2 = time.time()

		hypercube = np.array(hypercube)
		print('>>>>>>>>>>>>>')
		print("hypercube shape: ", hypercube.shape)
		print('>>>>>>>>>>>>>')

		unc = calunc(hypercube, [NAXIS1, NAXIS2, NAXIS3], \
			file_all+'_unc', hdr, wavALL)

		t_cal_unc = time.time()
		print("****** cal_unc_time = {:.0f} seconds ******".format(t_cal_unc - t2))
	else:
		unc = read_fits(file_all+'_unc')[1]
	
	fclean(path_tmp+'*.fits')


else:
	b2 = input("(Re)calculate uncertainty? [y/n] ")
	
	## Time offset
	t1 = time.time()

	cube0, wavALL = read_fits(file_all+'_0')[1:3]
	
	if b2=='y':
		hypercube=[]
		Nmc += 1
		for j in range(Nmc):
			if j==0:
			else:
				cube = read_fits(file_all+'_'+str(j))[1]
				hypercube.append(cube)
		
		hypercube = np.array(hypercube)
		print('>>>>>>>>>>>>>')
		print("hypercube shape: ", hypercube.shape)
		print('>>>>>>>>>>>>>')
		unc = calunc(hypercube, [NAXIS1, NAXIS2, NAXIS3], \
			file_all+'_unc', hdr, wavALL)

	else:
		unc = read_fits(file_all+'_unc')[1]

print('\n>>>>>>>>>>>>>>>>>>\n')
print("Final cube shape: ", unc.shape)
print('\n>>>>>>>>>>>>>>>>>>\n')
# for x in range(NAXIS1):
# 	for y in range(NAXIS2):
# 		if all(np.isnan(a)==0 for a in cube0[:,y,x]):
# 			specview(wvl0, [cube0[:,y,x]], [unc[:,y,x]], \
# 				savename=fig_path+'({}, {}).png'.format(x, y))

t_total = time.time()
print("****** total_time = {:.0f} seconds ******".format(t_total - t1))
