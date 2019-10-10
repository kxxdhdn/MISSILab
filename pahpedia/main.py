#!/usr/bin/env python
# -*- coding: utf-8 -*-

import time
t0 = time.time()

from astropy import units as u
import numpy as np
import matplotlib.pyplot as plt
from utils.myfunc import celest2deg, fclean
from utils.rwfits import *
from utils.impro import *
from utils.mapping import *
from utils.mc import calunc
from smooth import crop3D, choker, do_conv

source = 'N44'
chnl = ['SL2', 'SL1', 'LL2', 'LL1']

data_path = 'data/'+source+'/'
tmp_path = data_path+'tmp/'
out_path = data_path+'output/'
fig_path = 'figures/'+source+'/'
ref_file = data_path+source+'_LL1'

N_mc = 3

## crop config
# ra, dec = celest2deg(13., 37., 1., -29., 51., 55.5) # M83
# ra, dec = celest2deg(9., 55., 45.5, 69., 42., 16.2) # M82_N
# ra, dec = celest2deg(9., 55., 59.1, 69., 39., 21.8) # M82_S
# ra, dec = celest2deg(9., 55., 38.5, 69., 43., 58.4) # M82_O
ra, dec = celest2deg(9., 55., 38.5, 69., 43., 58.4) # N44
dx, dy = 15, 15

## read reprojection reference
# hdREF = read_fits(ref_file, wvl_mod=1)[2] # 3D nocrop
# w, hdREF = WCSextract(ref_file) # 2D nocrop
## slice & crop cube
crop3D(ref_file, tmp_path+'ref', \
	None, (ra, dec), (dx, dy), 1)[0]
w, hdREF = WCSextract(tmp_path+'ref_0000') # 2D cropped header
hdREF['EQUINOX'] = 2000.0
NAXIS1 = hdREF['NAXIS1']
NAXIS2 = hdREF['NAXIS2']

t_init = time.time()
print("****** init_time = {:.0f} seconds ******".format(t_init - t0))

b0 = input("(Re)do homegeneisation? [y/n] ")
if b0=='y':
	b1 = input("(Re)do Monte-Carlo? [y/n] ")
	if b1=='y':
		hypercube=[]
		N_mc += 1
	else:
		N_mc=1
	
	t1 = time.time()

	for j in range(N_mc):

		for i, ch in enumerate(chnl):

			data_filename = source+'_'+ch
			unc_filename = data_filename+'_unc'
			back_filename = source+'_off_'+ch
			filIN = data_path+data_filename
			filOUT = tmp_path+data_filename
			
			## slice cube
			if j==0:
				wvl = cubislice(filIN, filOUT, \
					None, None, 1)[0]
					# None, data_path+back_filename, 1)[0]
			else:
				wvl = cubislice(filIN, filOUT, \
					data_path+unc_filename, None, 1)[0]
					# data_path+unc_filename, data_path+back_filename, 1)[0]
			SList = []
			for k in range(np.size(wvl)):
				SList.append(filOUT+'_'+'0'*(4-len(str(k)))+str(k))
			
			## smooth
			choker(SList, wvl)
			do_conv()

			## reproject convolved cube
			cube = []
			for SLout in SList:
				data = read_fits(SLout+'_conv', False)[0]
				reproj, ft = rpj(SLout+'_conv', \
					SLout+'reproj', hdREF, write=0)
				cube.append(reproj)
			cube = np.array(cube)

			## extend wvl
			if i==0:
				wvl0 = wvl
				cubi = cube
			else:
				for k, lam in enumerate(wvl):
					if wvl0[-1]>lam:
						wvl0 = wvl0[:-1]
						cubi = cubi[:-1,:,:]
					else:
						break
				wvl0 = np.concatenate((wvl0, wvl[k+1:]))
				cubi = np.concatenate((cubi, cube[k+1:,:,:]))
				print("{} edge wvl deleted.".format(k*2))
			print(wvl0.shape)
			print(cubi.shape)
		NAXIS3 = np.size(wvl0)

		write_fits(out_path+source+'_'+str(j), cubi, hdREF, wvl0)
		if j==0:
			## no uncertainty added cube
			cube0 = cubi

			t_no_mc = time.time()
			print("****** no_mc_time = {:.0f} seconds ******".format(t_no_mc - t1))
		else:
			## ieme iteration finished
			print("----------------{}----------------".format(j))

			## append hypercube
			hypercube.append(cubi)
	if b1=='y':
		hypercube = np.array(hypercube)
		print("hypercube shape: ", hypercube.shape)
			
		## calculate MC uncertainty
		t2 = time.time()

		unc = calunc(hypercube, [NAXIS1, NAXIS2, NAXIS3], \
			out_path+source+'_unc', hdREF, wvl0)
		print("unc shape: ", unc.shape)

		t_cal_unc = time.time()
		print("****** cal_unc_time = {:.0f} seconds ******".format(t_cal_unc - t2))
	else:
		unc = read_fits(out_path+source+'_unc')[0]
	
	fclean(tmp_path+'*.fits')

else:
	b2 = input("(Re)calculate uncertainty? [y/n] ")
	if b2=='y':
		hypercube=[]
		N_mc += 1
		for j in range(N_mc):
			if j==0:
				cube0, wvl0 = read_fits(out_path+source+'_0')[0:2]
				NAXIS3 = np.size(wvl0)
			else:
				cubi = read_fits(out_path+source+'_'+str(j))[0]
				hypercube.append(cubi)
		hypercube = np.array(hypercube)
		print("hypercube shape: ", hypercube.shape)

		unc = calunc(hypercube, [NAXIS1, NAXIS2, NAXIS3], \
			out_path+source+'_unc', hdREF, wvl0)
		print("unc shape: ", unc.shape)
	else:
		cube0, wvl0 = read_fits(out_path+source+'_0')[0:2]
		unc = read_fits(out_path+source+'_unc')[0]

for x in range(NAXIS1):
	for y in range(NAXIS2):
		if all(np.isnan(a)==0 for a in cube0[:,y,x]):
			specview(wvl0, [cube0[:,y,x]], [unc[:,y,x]], \
				savename=fig_path+'({}, {}).png'.format(x, y))

t_total = time.time()
print("****** total_time = {:.0f} seconds ******".format(t_total - t_init))

# plt.show()
