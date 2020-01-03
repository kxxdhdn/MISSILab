#!/usr/bin/env python
# -*- coding: utf-8 -*-

import time
t0 = time.time()

import os
import numpy as np
import matplotlib.pyplot as plt
## astylo
from astylo.sinout import read_fits, write_fits, WCSextract, read_ascii
from astylo.processim import islice, icrop, iconvolve, iproject, wclean
from astylo.calib import intercalib, phot2phot
from astylo.myfunclib import fclean, MCerror, pix2sr
from astylo.splot import plot2d

##---------------------------
##       Initialisation
##---------------------------

Nmc = 2

src = 'M82'
## Not useful if need to modify IDL/conv_prog.pro
# src = input("Input source name: ")

instr = ['IRS', 'IRC']
# chnl = ['SL2', 'SL1', 'LL2'] # hdr['APERNAME'] related
chnl = ['SL2', 'SL1', 'LL2', 'LL1']

verbose = False

##---------------------------
##         Path Store
##---------------------------

## Data dir
path_data = '/Users/dhu/Data/PAHpedia/'+src+'/' ###
path_tmp = '/Users/dhu/Data/PAHpedia/tmp/' ### See also IDL/convolve_image.pro
if not os.path.exists(path_tmp):
	os.makedirs(path_tmp)
path_test = '/Users/dhu/Data/PAHpedia/tests/' ###
path_out = path_data+'output/' ###vvv See also
if not os.path.exists(path_out):
	os.makedirs(path_out)
file_ref = path_data+src+'_'+chnl[-1] # => project_ref
project_ref = path_out+src+'_ref' ###
file_all = path_out+src

## Current dir
path_cur = os.getcwd()+'/'
path_idl = path_cur+'/IDL/'

## Kernel list stocked in csv file
path_ker = '/Users/dhu/Data/kernels/' ###
cker = path_out+src+'_kernels' ### See also IDL/conv_prog.pro

## Intercalibration
path_phot = '/Users/dhu/Data/photometry/'+src+'/' ###
path_calib = path_data+'calib/' ###
if not os.path.exists(path_calib):
	os.makedirs(path_calib)
phot_s = 'IRAC4' ### phot for calib spec [spec2phot]
file_phot_s = path_phot+src+'_'+phot_s # => file_calib_s
Uconvert_s = True
file_calib_s = path_calib+src+'_'+phot_s
calib_s = path_calib+src+'_IRS_to_'+phot_s
phot_p = 'IRAC4_SINGS' ### phot for calib phot [phot2phot]
file_phot_p = path_phot+src+'_'+phot_p # => file_calib_p
Uconvert_p = False
file_calib_p = path_calib+src+'_'+phot_p

##---------------------------
##         Read Data
##---------------------------

## Image
##-------
coord = read_ascii(path_cur+'data/coord')
for c in coord:
	if c[0]==src:
		ra, dec = float(c[1]), float(c[2])
dx, dy = 15, 15 ###

## Kernel
##--------
psf = [2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6.]
kernelist = []
psf_ref = 'MIPS1' ###

if psf_ref=='Gauss_6.0':
	file_ker = '_to_Gauss_06.0'
elif psf_ref=='IRAC3': # 2.11 (no LL1)
	file_ker = '_IRAC_5.8'
elif psf_ref=='IRAC4': # 2.82 (no LL1)
	file_ker = '_IRAC_8.0'
elif psf_ref=='MIPS1': # 6.43
	file_ker = '_to_MIPS_24'
elif psf_ref=='WISE3': # 6.60
	file_ker = '_WISE_MAP_11.6'
elif psf_ref=='WISE4': # 11.89
	file_ker = '_WISE_MAP_22.1'
else:
	file_ker = None

if file_ker is not None:
	for p in psf:
		kernelist.append(path_ker+'Kernel_HiRes_Gauss_0'+str(p)+file_ker) ###

##---------------------------
##         Processing
##---------------------------

## Cropped ref image
##-------------------
icrop(file_ref, project_ref, cenval=(ra,dec), sizpix=(dx,dy), wmod=1)
hdr = read_fits(project_ref)[0]
NAXIS1 = hdr['NAXIS1']
NAXIS2 = hdr['NAXIS2']

t_init = time.time()
print(">> init_time = {:.0f} seconds <<".format(t_init - t0))

b0 = input("(Re)do homegeneisation? [y/n] ")
# b0 = 'y'
if b0=='y':
	b1 = input("(Re)do Monte-Carlo? [y/n] ")
	
	## Time offset
	t1 = time.time()

	# b1 = 'y'
	if b1=='y':
		hypercube=[]
	else:
		Nmc = 0
	
	for j in range(Nmc+1):
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
			# 	sl = islice(file_data, file_slice)
			# 	wvl = sl.wave()
			# 	slist = sl.filenames()
			# else: # add unc
			# 	islice(filIN=file_data, filSL=file_slice, uncIN=file_unc)

			if file_ker is not None:
				## Smooth images [IDL]
				##---------------------
				if j==0: # unc not added
					conv = iconvolve(filIN=file_data, filKER=kernelist, \
						saveKER=cker, wmod=1, \
						filTMP=file_slice, filOUT=file_conv)
				else: # add unc
					conv = iconvolve(filIN=file_data, filKER=kernelist, \
						saveKER=cker, wmod=1, uncIN=file_unc, \
						filTMP=file_slice, filOUT=file_conv)

				conv.do_conv(ipath=path_idl)

				wavALL.extend(conv.wave())

				file_proj = file_conv
			else:
				file_proj = file_data

			## Reproject convolved cube
			##--------------------------
			pro = iproject(filIN=file_proj, filREF=project_ref, \
				filTMP=file_slice)#, filOUT=file_out+'_'+str(j))
			
			slices.extend(pro.filenames())
			cubi.append(pro.image()) # Here, pro.image is a 3D cube
		
			fclean(file_conv+'.fits')

		cube = np.concatenate(cubi)

		# hdr['APERNAME'] = 'SL2+SL1+LL2(ref)'
		hdr['APERNAME'] = 'SL2+SL1+LL2+LL1(ref)'
		comment = "Homegeneized cube produced by [SPEXTRACT] routine. "
		write_fits(file_all+'_'+str(j), hdr, cube, wavALL, 1, COMMENT=comment)

		if j==0:
			## no uncertainty added cube
			cube0 = cube

			t_no_mc = time.time()
			print(">> no_mc_time = {:.0f} seconds <<".format(t_no_mc - t1))
		else:
			## append hypercube
			hypercube.append(cube)

			## ieme iteration finished
			print("---------------- {} ----------------".format(j))
	
	if verbose==False:
		fclean(path_tmp+'*.fits')
	b2 = 'y'
else:
	b2 = input("(Re)calculate uncertainty? [y/n] ")

## Calculate MC uncertainties
##----------------------------
## Time offset
t2 = time.time()

cube0, wavALL = read_fits(file_all+'_0')[1:3]

if b2=='y':
	hypercube=[]
	for j in range(Nmc+1):
		if j==0:
			pass
		else:
			cube = read_fits(file_all+'_'+str(j))[1]
			hypercube.append(cube)

	hypercube = np.array(hypercube)
	print('>>>>>>>>>>>>>')
	print("hypercube shape: ", hypercube.shape)
	print('>>>>>>>>>>>>>')

	unc = MCerror(hypercube)
	write_fits(file_all+'_unc', hdr, unc, wavALL, 1)

	t_cal_unc = time.time()
	print(">> cal_unc_time = {:.0f} seconds <<".format(t_cal_unc - t2))
else:
	unc = read_fits(file_all+'_unc')[1]

## Clean Wavelengths
##-------------------
cube0, wavALL = wclean(file_all+'_0', wmod=1, filOUT=file_all, verbose=True)
unc = wclean(file_all+'_unc', wmod=1, filOUT=file_all+'_unc')[0]

print('\n>>>>>>>>>>>>>>>>>>\n')
print("Final cube shape: ", unc.shape)
print('\n>>>>>>>>>>>>>>>>>>\n')

##---------------------------
##      Intercalibration
##---------------------------

## Convert Jy/pix to MJy/sr (Optional)
##-------------------------------------
hdr_s, im_s = read_fits(file_phot_s)[0:2]
# if hdr_s['SIGUNIT']=='Jy/pix              / Unit of the map': # DustPedia
if Uconvert_s==True:
	im_s = im_s * 1.e-6 / pix2sr(1., hdr_s['CDELT1'])
	hdr_s['SIGUNIT'] = 'MJy/sr'
write_fits(file_calib_s, hdr_s, im_s)

hdr_p, im_p = read_fits(file_phot_p)[0:2]
# if hdr_p['SIGUNIT']=='Jy/pix              / Unit of the map': # DustPedia
if Uconvert_p==True:
	im_p = im_p * 1.e-6 / pix2sr(1., hdr_p['CDELT1'])
	hdr_p['SIGUNIT'] = 'MJy/sr'
write_fits(file_calib_p, hdr_p, im_p)

## Crop to save time
##-------------------
# dx_phot = int(2. * dx)
# dy_phot = int(2. * dy)
# cro_s = icrop(file_calib_s, file_calib_s, \
# 	cenval=(ra, dec), sizpix=(dx_phot, dy_phot))
# cro_p = icrop(file_calib_p, file_calib_p, \
# 	cenval=(ra, dec), sizpix=(dx_phot, dy_phot))

## spec2phot
##-----------

## phot2phot
##-----------
## Use reprojection to crop
pro_s = iproject(file_calib_s, file_all, filOUT=file_calib_s)
p2p = phot2phot(filIN=file_calib_p, \
	filREF=file_calib_s, filOUT=file_calib_p)
newim_p = p2p.image()
newim_s = pro_s.image()

plt.figure()
x = np.arange(1.,1.e4,1.)
plt.plot(x, x, c='k')
plt.scatter(newim_p, newim_s, c='m', s=1., label='p2p')
plt.xscale('symlog')
plt.yscale('symlog')
# plt.xlim((0,2.e2))
# plt.ylim((0,2.e2))
plt.xlabel(phot_p)
plt.ylabel(phot_s)
plt.legend()

##---------------------------
##           plot
##---------------------------

# plot2d(wavALL, cube0[:,10,10], yerr=unc[:,10,10])

plt.show()

t_total = time.time()
if b0=='y':
	print(">> total_time = {:.0f} seconds <<".format(t_total - t1))
else:
	print(">> total_time = {:.0f} seconds <<".format(t_total - t2))
