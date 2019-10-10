#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

matplotlib applications

"""
from astropy import units as u
import numpy as np
import matplotlib.pyplot as plt
from utils.myfunc import pix2sr, fclean
from utils.rwfits import *
from utils.synthetic_photometry import *
from utils.impro import *
from utils.mapping import *
from utils.mc import calunc
from smooth import marker, choker, do_conv

data_path = 'data/N44/output/'
data_filename = 'N44_0'
unc_filename = 'N44_unc'
# back_filename = 'm83_off_LL1'
out_path = data_path+'../convolved/'
rpj_path = data_path+'../reprojection/'

ph_path = '/Users/dhu/data/photometry/LMC/'
ph_filename = 'SAGEA_LMC_8x8_mosaic_I3'
# ph_filename = 'NGC5236_Spitzer_5.8'
# ph_filename = 'NGC5236_Spitzer_8.0'
# ph_filename = 'NGC3034_Spitzer_5.8'
# ph_filename = 'NGC3034_Spitzer_8.0'
out_ph = data_filename+'_'+ph_filename

N_mc = 3

## uncertainty adjustment factor for SL2, SL1, LL2, LL1 (IRS Handbook P.175)
IRSunc = [.011, .011, .016, .15]
## calibration uncertainty for 24, 70, 160 microns (MIPS Handbook P.79)
MIPSunc = [.04, .07, .12]
## calibration uncertainty for 3.6, 4.5, 5.8, 8.0 microns (IRAC Handbook P.36)
IRACunc = [.03, .03, .03, .03]

## 3D cube slicing
wvl, hdr = cubislice(data_path+data_filename, out_path+data_filename, \
	None, None, 0)
	# None, data_path+back_filename, '_')
NAXIS1 = hdr['NAXIS1']
NAXIS2 = hdr['NAXIS2']
NAXIS3 = np.size(wvl)
SList = []
for k in range(NAXIS3):
	sl_file = out_path+data_filename+'_'+'0'*(4-len(str(k)))+str(k)
	SList.append(sl_file)

## smooth
# choker(SList, wvl)
# do_conv()

## rebuild convolved cube
cube = []
for SLout in SList:
	# cubi = read_fits(SLout+'_conv', False)[0]
	cubi = read_fits(SLout, False)[0]
	cube.append(cubi)
cube = np.array(cube)

## photometry data
dx, dy = NAXIS1*25, NAXIS2*25
ra, dec = hdr['CRVAL1'], hdr['CRVAL2']
crop(ph_path+ph_filename, data_path+out_ph, \
	(ra, dec), (dx, dy))
crop(ph_path+ph_filename+'_Error', data_path+out_ph+'_unc', \
	(ra, dec), (dx, dy))
## [Alternative] hextract crop
# hextract(ph_path+ph_filename, data_path+out_ph, \
	# 9250, 10750, 5500, 7000)
# mSubimage(infile=ph_path+ph_filename+'.fits', outfile=data_path+out_ph+'.fits', ra=ra, dec=dec, xsize=dx, ysize=dy, mode=1)

## convert Jy/pix to MJy/sr
old_im, PhDr = read_fits(data_path+out_ph, False)
new_im = old_im * 1e-6 / pix2sr(1., PhDr['CDELT1'])
write_fits(data_path+out_ph, new_im, PhDr)
old_im, PhDr = read_fits(data_path+out_ph+'_unc', False)
new_im = old_im * 1e-6 / pix2sr(1., PhDr['CDELT1'])
write_fits(data_path+out_ph+'_unc', new_im, PhDr)
## convolve photometry file
marker(data_path+out_ph)
do_conv()
## mips024 -> irsll1 reprojection
photometry, ft = rpj(data_path+out_ph+'_conv', rpj_path+out_ph, hdr)
PHerr, ft = rpj(data_path+out_ph+'_unc_conv', rpj_path+out_ph, hdr)

w = WCSextract(data_path+data_filename)[0]

## inter-calib
Fnu_filt = []
Unc_filt = []
Fnu_ph = []
Unc_ph = []
abpts = []
## do synthetic photometry
cube_m = np.zeros((NAXIS3+1)*NAXIS2*NAXIS1).reshape((NAXIS3+1), NAXIS2, NAXIS1)
## correction of the interpolation (at the beginning of spectra)
wvl_m = [wvl[0]-3.]
wvl_m.extend(wvl)
wvl_m = np.array(wvl_m)
cube_m[1:,:,:] = np.copy(cube)
# wcen, data0, sig = synthetic_photometry(wvl_m, cube_m, ['MIPS1'])
wcen, data0, sig = synthetic_photometry(wvl_m, cube_m, ['IRAC4'])

## save data
write_fits(rpj_path+data_filename+'_calib', data0, hdr)

"""
----------------------------------------
Uncertainty propagation by Monte-Carlo
----------------------------------------

"""
data_mc = []
b0 = input("(Re)do Monte-Carlo? [y/n] ")
for i in range(N_mc):
	if b0=='y':
		## 3D cube slicing with uncertainty added
		cubislice(data_path+data_filename, out_path+data_filename, \
			data_path+unc_filename, None, 0)
			# data_path+unc_filename, data_path+back_filename, '_')

		## smooth
		# choker(SList, wvl)
		# do_conv()
		# fclean(data_path+'../convolved/*_.fits')

		## rebuild convolved cube
		cube = []
		for SLout in SList:
			# cubi = read_fits(SLout+'conv', False)[0]
			cubi = read_fits(SLout, False)[0]
			cube.append(cubi)
		cube = np.array(cube)

		## inter-calib
		cube_m = np.zeros((NAXIS3+1)*NAXIS2*NAXIS1).reshape((NAXIS3+1), NAXIS2, NAXIS1)
		## do synthetic photometry
		cube_m[1:,:,:] = np.copy(cube)
		# data = synthetic_photometry(wvl_m, cube_m, ['MIPS1'])[1]
		data = synthetic_photometry(wvl_m, cube_m, ['IRAC4'])[1]
		## save data
		write_fits(rpj_path+data_filename+'_calib_'+'0'*(4-len(str(i)))+str(i), data, hdr)
		
		## ieme iteration finished
		print("----------------{}----------------".format(i+1))
	else:
		data = read_fits(rpj_path+data_filename+'_calib_'+'0'*(4-len(str(i)))+str(i), False)[0]
	data_mc.append(data)
data_mc = np.array(data_mc)
FILTerr = calunc(data_mc, [NAXIS1, NAXIS2])

fclean(data_path+'../convolved/*_conv.fits')

# reprojection mask & correlation linearity check
for x in range(NAXIS1):
	for y in range(NAXIS2):
		if data0[y,x]!=0 and photometry[y,x]!=np.nan:
			Fnu_filt.append(data0[y,x])
			Unc_filt.append(np.sqrt(FILTerr[y,x]**2 + (data0[y,x]*IRSunc[3])**2))
			Fnu_ph.append(photometry[y,x])
			# Unc_ph.append(np.sqrt(PHerr[y,x]**2 + (photometry[y,x]*MIPSunc[0])**2))
			Unc_ph.append(np.sqrt(PHerr[y,x]**2 + (photometry[y,x]*IRACunc[0])**2))
			if data0[y,x]/photometry[y,x]>2. or photometry[y,x]/data0[y,x]>2.:
				abpts.append([x,y])
		## mask data according to reprojection footprint
		else:
			data0[y,x]=np.nan

## print pixels that are beyond linearity
# print("aberrant pixels: ", abpts)
print("aberrant pix nb: ", np.array(abpts).shape)
y, x = 8, 8
specview(wvl_m, [cube_m[:,y,x]], None, \
	wfilt=wcen, Fnu=[data0[y,x], photometry[y,x]], filt_width=6.4)
factor = calibview(Fnu_ph, Fnu_filt, yerr=Unc_filt, xerr=Unc_ph, \
	threshold=2., \
	xlabel="IRAC4 (MJy/sr)", ylabel="IRS SL (MJy/sr)", \
	savename='figures/calib/'+data_filename+'_calib.png')
imview([data0, photometry], (1,2), w, add_pts=abpts, figsize=(12,6), rot=1, \
	savename='figures/calib/'+data_filename+'_field.png')

b1 = input("Correct spectra? [y/n] ")
if b1=='y':
	im = specorrect(data_path+data_filename, data_path+data_filename, \
		factor, 4., 14.5, wmod=0)

plt.show()
