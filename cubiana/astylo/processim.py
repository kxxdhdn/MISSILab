#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

PROCESS IMage

"""

import math
import numpy as np
from scipy.io import readsav
from reproject import reproject_interp
import subprocess as SP

## astylo
from sinout import read_fits, write_fits, WCSextract, read_csv, write_csv
from myfunclib import fclean, closest, bsplinterpol


def wmask(filIN, filOUT=None):
	'''
	MASK Wavelengths

	--- INPUT ---
	filIN       input fits file 
	filOUT      overwrite fits file (Default: NO)
	--- OUTPUT ---
	data_new    new fits data
	wave_new    new fits wave
	'''
	pass


def wclean(filIN, wmod=0, cmod='eq', cfile=None, \
	filOUT=None, verbose=False):
	'''
	CLEAN Wavelengths

	--- INPUT ---
	filIN       input fits file
	wmod        wave mode
	cmod        clean mode (Default: 'eq')
	cfile       input csv file (archived info)
	filOUT      overwrite fits file (Default: NO)
	verbose     display wclean info (Default: False)
	--- OUTPUT ---
	data_new    new fits data
	wave_new    new fits wave
	'''
	hdr, data, wave = read_fits(filIN, wmod=wmod)
	Nw = len(wave)
	
	ind = [] # list of indices of wvl to remove
	if cfile is not None:
		indarxiv = read_csv(cfile, 'Ind')[0]
		ind = []
		for i in indarxiv:
			ind.append(int(i))
	else:
		## Detect crossing wvl
		##---------------------
		for i in range(Nw-1):
			if wave[i]>=wave[i+1]: # found wave(i+1), i_max=Nw-2
				
				wmin = -1 # lower limit: closest wave smaller than wave[i+1]
				wmax = 0 # upper limit: closest wave larger than wave[i]
				
				for j in range(i+1):
					dw = wave[i+1] - wave[i-j]
					if dw>0: # found the closest smaller wave[i-j]
						wmin = i-j
						break # only the innermost loop
				if wmin==-1:
					print('WARNING: Left side fully covered! ')
				
				for j in range(Nw-i-1):
					dw = wave[i+1+j] - wave[i]
					if dw>0: # found the closest larger wave[i+1+j]
						wmax = i+1+j
						break
				if wmax==0:
					print('WARNING: right side fully covered! ')

				Nw_seg = wmax-wmin-1 # number of crossing wvl in segment
				wave_seg = [] # a segment (every detect) of wave
				ind_seg = [] # corresponing segment for sort use
				for k in range(Nw_seg):
					wave_seg.append(wave[wmin+1+k])
					ind_seg.append(wmin+1+k)
				## index list of sorted wave_seg
				ilist = sorted(range(len(wave_seg)), key=wave_seg.__getitem__)
				## index of wave_seg center
				icen = math.floor((Nw_seg-1)/2)

				## Visualisation (for test use)
				##------------------------------
				# print('wave, i: ', wave[i], i)
				# print('wave_seg: ', wave_seg)
				# print('ind_seg: ', ind_seg)
				# print('ilist: ', ilist)
				# print('icen: ', icen)

				## Remove all crossing wvl between two channels
				##----------------------------------------------
				if cmod=='all': # most conservative but risk having holes
					pass
				## Remove (almost) equal wvl (NOT nb of wvl!) for both sides
				##-----------------------------------------------------------
				elif cmod=='eq': # (default)
					## Select ascendant pair closest to segment center
					for k in range(icen):
						if ilist[icen]>ilist[0]: # large center
							if ilist[icen-k]<ilist[0]:
								for p in range(ilist[icen-k]+1):
									del ind_seg[0]
								for q in range(Nw_seg-ilist[icen]):
									del ind_seg[-1]
								break
						else: # small center
							if ilist[icen+k]>ilist[0]:
								for p in range(ilist[icen]+1):
									del ind_seg[0]
								for q in range(Nw_seg-ilist[icen+k]):
									del ind_seg[-1]
								break
				## Leave 2 closest wvl not crossing
				##----------------------------------
				elif cmod=='closest_left':
					for k in range(ilist[0]):
						del ind_seg[0]
				elif cmod=='closest_right':
					for k in range(Nw_seg-ilist[0]):
						del ind_seg[-1]
				## Others
				##--------
				else:
					print('ERROR: Not supported clean mode! ')

				# print('ind_seg (final): ', ind_seg)
				ind.extend(ind_seg)

	## Do clean
	##----------
	data_new = np.delete(data, ind, axis=0)
	wave_new = list(np.delete(np.array(wave), ind))

	## Display clean detail
	##----------------------
	if verbose==True:
		print('\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
		print('Number of wavelengths deleted: ', len(ind))
		print('Ind, wavelengths: ')
		for i in ind:
			print(i, wave[i])
		print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n')

	## Overwrite fits file
	##---------------------
	if filOUT is not None:
		# comment = 'Wavelength removal info in _wclean_info.csv'

		write_fits(filOUT, hdr, data_new, wave_new) # hdr auto changed
		
		## Write csv file
		wlist = []
		for i in ind:
			wlist.append([i, wave[i]])
		write_csv(filOUT+'_wclean_info', \
			header=['Ind', 'Wavelengths'], dataset=wlist)

	return data_new, wave_new

def interfill(arr, axis):
	'''
	FILL undersampling/artificial gap by (bspl)INTERpolation

	--- INPUT ---
	arr         array
	axis        axis along which interpolation
	--- OUTPUT ---
	newarr      new array
	'''
	print(">> fill gaps with b-splines <<")

	axsh = arr.shape
	NAXIS = np.size(axsh)
	newarr = np.copy(arr)
	if NAXIS==1: # 1D array
		x = np.arange(axsh[0])
		for i in range(axsh[0]):
			newarr = bsplinterpol(x, arr, x)
	if NAXIS==2: # no wavelength
		if axis==0: # col direction
			y = np.arange(axsh[0])
			for i in range(axsh[1]):
				col = bsplinterpol(y, arr[:,i], y)
				for j in range(axsh[0]):
					newarr[j,i] = col[j]
		elif axis==1: # row direction
			x = np.arange(axsh[1])
			for j in range(axsh[0]):
				row = bsplinterpol(x, arr[j,:], x)
				for i in range(axsh[1]):
					newarr[j,i] = row[i]
		else:
			print('ERROR: Unkown axis! ')
	elif NAXIS==3:
		if axis==0: # fill wavelength
			z = np.arange(axsh[0])
			for i in range(axsh[2]):
				for j in range(axsh[1]):
					wvl = bsplinterpol(z, arr[:,j,i], z)
					for k in range(axsh[0]):
						newarr[k,j,i] = wvl[k]
		elif axis==1: # col direction
			y = np.arange(axsh[1])
			for k in range(axsh[0]):
				for i in range(axsh[2]):
					col = bsplinterpol(y, arr[k,:,i], y)
					for j in range(axsh[1]):
						newarr[k,j,i] = col[j]
		elif axis==2: # row direction
			x = np.arange(axsh[2])
			for k in range(axsh[0]):
				for j in range(axsh[1]):
					row = bsplinterpol(x, arr[k,j,:], x)
					for i in range(axsh[2]):
						newarr[k,j,i] = row[i]
		else:
			print('ERROR: Unkown axis! ')
	else:
		print('ERROR: array shape not supported! ')

	return newarr

def hextract(filIN, filOUT, x0, x1, y0, y1):
	'''
	Crop 2D image with pixel sequence numbers
	[ref]
	IDL lib hextract
	https://idlastro.gsfc.nasa.gov/ftp/pro/astrom/hextract.pro
	'''
	oldim = read_fits(filIN)[1]
	hdr, w = WCSextract(filIN)[0:2]
	# hdr['NAXIS1'] = x1 - x0 + 1
	# hdr['NAXIS2'] = y1 - y0 + 1
	hdr['CRPIX1'] += -x0
	hdr['CRPIX2'] += -y0
	newim = oldim[y0:y1+1, x0:x1+1]

	write_fits(filOUT, hdr, newim)

	return newim

##-----------------------------------------------

##			"improve" based tools

##-----------------------------------------------

class improve:
	'''
	IMage PROcessing VEssel
	'''
	def __init__(self, filIN, wmod=0):
		
		## INPUTS
		self.filIN = filIN
		self.wmod = wmod

		## read image/cube
		## self.hdr is a 2D (reduced) header
		self.hdr, self.w, self.is3d = WCSextract(filIN)
		self.NAXIS1 = self.hdr['NAXIS1']
		self.NAXIS2 = self.hdr['NAXIS2']
		print("Raw size (pix): {} * {}".format(
			self.NAXIS1, self.NAXIS2))
		## 3D cube slicing
		if self.is3d==True:
			self.im, self.wvl = read_fits(filIN, wmod=wmod)[1:3]
			self.NAXIS3 = len(self.wvl)
		else:
			self.im = read_fits(filIN, wmod=wmod)[1]
			self.wvl = None

	def addunc(self, uncIN=None, wmod=0):
		## uncertainty adding
		mu, sigma = 0., 1.

		if uncIN is not None:
			unc = read_fits(uncIN, wmod)[1]
			theta = np.random.normal(mu, sigma, unc.shape)
			## unc has the same dimension with im
			self.im += theta * unc

	def slice(self, filSL, suffix=None):
		## 3D cube slicing
		if self.is3d==True:
			slicnames = []
			for k in range(self.NAXIS3):
				## output filename list
				f = filSL+'_'+'0'*(4-len(str(k)))+str(k)
				if suffix is not None:
					f += suffix
				slicnames.append(f)
				comment = "NO.{} image [SLICE]d from {}.fits".format(k, self.filIN)
				write_fits(f, self.hdr, self.im[k,:,:])
		else:
			print('Input file is a 2D image which cannot be sliced! ')
			exit()

		return slicnames
	
	def rectangle(self, CENval, SIZpix, CENpix=None, filOUT=None):
		## rectangle center in (ra, dec) if CENpix is None
		## rectangle size in pix
		## [optional] 2D raw data to crop
		
		if CENpix is None:
			## convert coord
			xcen, ycen = self.w.all_world2pix(CENval[0], CENval[1], 1)
			print("Crop centre (ra, dec): [{:.8}, {:.8}]".format(*CENval))
			if not (0<xcen<self.NAXIS1 and 0<ycen<self.NAXIS2):
				print("ERROR: crop centre overpassed image border! ")
				exit()
		else:
			xcen, ycen = CENpix[0], CENpix[1]
			print("Crop centre (x, y): [{}, {}]".format(*CENpix))
		
		print("Crop size (pix): [{}, {}].".format(*SIZpix))
		
		## lowerleft origin
		xmin = math.floor(xcen - SIZpix[0]/2.)
		ymin = math.floor(ycen - SIZpix[1]/2.)
		xmax = xmin + SIZpix[0]
		ymax = ymin + SIZpix[1]
		if not (xmin>=0 and xmax<=self.NAXIS1 and ymin>=0 and ymax<=self.NAXIS2):
			print("ERROR: crop region overpassed image border! ")
			exit()

		## OUTPUTS
		##---------
		## new image
		if self.is3d==True:
			self.im = self.im[:, ymin:ymax, xmin:xmax]
			## recover 3D non-reduced header
			self.hdr = read_fits(self.filIN, self.wmod)[0]
		else:
			self.im = self.im[ymin:ymax, xmin:xmax]
		## modify header
		self.hdr['CRPIX1'] += -xmin
		self.hdr['CRPIX2'] += -ymin
		## write cropped image/cube
		if filOUT is not None:
			comment = "[RECTANGLE] cropped at centre: [{:.8}, {:.8}]. ".format(*CENval)
			# comment = "with size [{}, {}] (pix).".format(*SIZpix)

			write_fits(filOUT, self.hdr, self.im, wave=self.wvl, \
				COMMENT=comment)

		return self.im

class slicube(improve):
	'''
	SLIce CUBE
	'''
	def __init__(self, filIN, filSL, wmod=0, \
		uncIN=None, wmod_unc=0, suffix=None):
		super().__init__(filIN, wmod)
		
		self.addunc(uncIN, wmod_unc)

		self.slicnames = self.slice(filSL, suffix)

	def image(self):
		return self.im

	def wave(self):
		return self.wvl

	def slice_names(self):
		return self.slicnames

class crop(improve):
	'''
	CROP 2D image or 3D cube
	'''
	def __init__(self, filIN, CENval, SIZpix, \
		wmod=0, uncIN=None, wmod_unc=0, \
		shape='box', filOUT=None):
		## slicrop: slice 
		super().__init__(filIN, wmod)

		self.addunc(uncIN, wmod_unc)
		
		if shape=='box':
			self.cropim = self.rectangle(
				CENval=CENval, SIZpix=SIZpix, filOUT=filOUT)

	def image(self):
		return self.cropim

	def wave(self):
		return self.wvl

class project(improve):
	'''
	PROJECT 2D image or 3D cube
	'''
	def __init__(self, filIN, filREF=None, hdREF=None, wmod=0, \
		uncIN=None, wmod_unc=0, filTMP=None, filOUT=None):
		super().__init__(filIN, wmod)
		## input hdREF must be (reduced) 2D header

		self.addunc(uncIN, wmod_unc)

		if filREF is not None:
			hdREF = WCSextract(filREF)[0]
			hdREF['EQUINOX'] = 2000.0
		
		if hdREF is not None:
			if self.is3d==True:
				if filTMP is not None:
					filSL = filTMP
				else:
					filSL = self.filIN
				slicIN = self.slice(filSL=filSL, suffix='_')
				im = []
				self.slicnames = []
				for f in slicIN:
					im.append(reproject_interp(f+'.fits', hdREF)[0])
					if filTMP is not None:
						write_fits(f+'_pro', hdREF, im)
						self.slicnames.append(f+'_pro')
					fclean(f+'.fits')

				self.projim = np.array(im)
				## recover non-reduced 3D header
				if filREF is not None:
					hdr = read_fits(filREF, self.wmod)[0]
			else:
				self.projim, footprint = reproject_interp(filIN+'.fits', hdREF)
		else:
			print('ERROR: Can not find projection reference! ')
			
		if filOUT is not None:
			comment = "Re[PROJECT]ed image. "
			
			if self.is3d==True:
				write_fits(filOUT, hdr, self.projim, wave=self.wvl, \
					COMMENT=comment)
			else:
				write_fits(filOUT, hdREF, self.projim, wave=self.wvl, \
					COMMENT=comment)
	
	def image(self):
		return self.projim

	def wave(self):
		return self.wvl

	def slice_names(self):
		return self.slicnames

class iconvolve(improve):
	'''
	(IDL based) CONVOLVE 2D image or 3D cube with given kernels

	--- INPUT ---
	filKER      convolution kernel(s) (tuple or list)
	--- OUTPUT ---
	'''
	def __init__(self, filIN, filKER, saveKER, \
		wmod=0, uncIN=None, wmod_unc=0, cfile=None, \
		psf=None, filTMP=None, filOUT=None):
		## INPUTS
		super().__init__(filIN, wmod)
		
		self.addunc(uncIN, wmod_unc)

		## input kernel file list
		self.filKER = filKER
		## doc (csv) file of kernel list
		self.saveKER = saveKER
		self.cfile = cfile
		self.filTMP = filTMP
		self.filOUT = filOUT

		## INIT
		if psf is None:
			self.psf = [2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.]
		else:
			self.psf = psf
		self.sigma_lam = None
				
	def spitzer_irs(self):
		'''
		Spitzer/IRS PSF profil
		[ref]
		Pereira-Santaella, Miguel, Almudena Alonso-Herrero, George H.
		Rieke, Luis Colina, Tanio Díaz-Santos, J.-D. T. Smith, Pablo G.
		Pérez-González, and Charles W. Engelbracht. “Local Luminous
		Infrared Galaxies. I. Spatially Resolved Observations with the
		Spitzer Infrared Spectrograph.” The Astrophysical Journal
		Supplement Series 188, no. 2 (June 1, 2010): 447.
		doi:10.1088/0067-0049/188/2/447.
		'''
		sim_par_wave = [0, 13.25, 40.]
		sim_par_fwhm = [2.8, 3.26, 10.1]
		sim_per_wave = [0, 15.5, 40.]
		sim_per_fwhm = [3.8, 3.8, 10.1]
		
		## fwhm (arcsec)
		fwhm_par = np.interp(self.wvl, sim_par_wave, sim_par_fwhm)
		fwhm_per = np.interp(self.wvl, sim_per_wave, sim_per_fwhm)
		#fwhm_lam = np.sqrt(fwhm_par * fwhm_per)
		
		## sigma (arcsec)
		sigma_par = fwhm_par / (2. * np.sqrt(2.*np.log(2.)))
		sigma_per = fwhm_per / (2. * np.sqrt(2.*np.log(2.)))
		self.sigma_lam = np.sqrt(sigma_par * sigma_per)
		
	def choker(self, filIM):
		## CHOose KERnel(s)

		if self.cfile is not None: # read archived info
			karxiv = read_csv(self.cfile, 'Images', 'Kernels')
			klist = []
			for k in karxiv:
				klist.append(k)
		else: # create list
			klist = []
			for i, image in enumerate(filIM):
				## check PSF profil (or is not a cube)
				if self.sigma_lam is not None:
					image = filIM[i]
					ind = closest(self.psf, self.sigma_lam[i])
					kernel = self.filKER[ind]
				else:
					image = filIM[0]
					kernel = self.filKER[0]
				## klist line elements: image, kernel
				k = [image, kernel]
				klist.append(k)

			## write csv file
			write_csv(self.saveKER, header=['Images', 'Kernels'], dataset=klist)

	def do_conv(self, ipath):
		
		if self.is3d==True:
			if self.filTMP is not None:
				filSL = self.filTMP
			else:
				filSL = self.filIN
			
			slicIN = self.slice(filSL)
			self.spitzer_irs()
			self.choker(slicIN)
		else:
			if self.filTMP is not None: # ?? TO DELETE
				filIN = [self.filTMP]
			else:
				filIN = [self.filIN]
			
			self.choker(filIN)

		SP.call('cd '+ipath+'\nidl conv.pro', shell=True)

		## OUTPUTS
		##---------
		if self.is3d==True:
			im = []
			self.slicnames = []
			for f in slicIN:
				im.append(read_fits(f+'_conv')[1])
				self.slicnames.append(f+'_conv')

			self.convim = np.array(im)
			## recover non-reduced 3D header
			self.hdr = read_fits(self.filIN, self.wmod)[0]
		else:
			self.convim = read_fits(self.filIN+'_conv')[1]
		
		if self.filOUT is not None:
			comment = "[ICONVOLV]ed by G. Aniano's IDL routine."
			# comment = 'https://www.astro.princeton.edu/~ganiano/Kernels.html'
			write_fits(self.filOUT, self.hdr, self.convim, wave=self.wvl, \
				COMMENT=comment)

	def image(self):
		return self.convim

	def wave(self):
		return self.wvl

	def slice_names(self):
		return self.slicnames

class slitextract(improve):
	'''
	AKARI/IRC spectroscopy slit coord extraction

	--- INPUT ---
	savIN       input sav file (IDL)
	filIRC      IRC N3 (long exp) frame corrected (2MASS corrected)
	--- OUTPUT ---
	'''
	def __init__(self, savIN, filIRC, filOUT):
		super().__init__(filIRC)

		im = readsav(savIN+'.sav', python_dict=True)['specimage_n_wc']
		im = im[::-1] # -> ascending order
		wave = readsav(savIN+'.sav', python_dict=True)['wave_array']
		wave = wave[::-1] # -> ascending order
		table = readsav(savIN+'.sav', python_dict=True)['source_table']
		ref_x = table['image_x'][0] # slit ref x
		ref_y = table['image_y'][0] # slit ref y
		dx = im.shape[1] # slit length
		dy = 1 # slit width: 256*5"/10' for 'Ns' & 256*3"/10' for 'Nh' 
		cube = im.reshape(im.shape[0], dy, dx)

		## Modify cube header
		self.rectangle(CENval=None, SIZpix=(dx, dy), CENpix=(ref_x, ref_y))
		self.hdr['CTYPE3'] = 'WAVELENGTH'
		self.hdr['CUNIT1'] = 'deg'
		self.hdr['CUNIT2'] = 'deg'
		self.hdr['BUNIT'] = 'MJy/sr'
		self.hdr['EQUINOX'] = 2000.0
		comment = "Assembled AKARI/IRC slit spectroscopy cube. "

		write_fits(filOUT, self.hdr, cube, wave=wave, \
				COMMENT=comment)

class imontage(improve):
	'''
	Image/cube mosaic

	--- INPUT ---
	--- OUTPUT ---
	'''
	def __init__(self, arg):
		super().__init__(filIN, wmod)

		self.arg = arg
		

"""
------------------------------ MAIN (test) ------------------------------
"""
if __name__ == "__main__":

	from myfunclib import MCerror

	import os
	path = os.getcwd() + '/../tests/'

	## Test slitextract
	##------------------
	savIN = path + 'data/1420415.1.N3_NG.IRC_SPECRED_OUT'
	filIRC = path + 'data/F011213824_N002'
	filOUT = path + 'out/M83_AKARI'

	slitextract(savIN, filIRC, filOUT)

	## Test interfill
	##----------------
	# Nmc = 2

	# filIN = path + 'data/IC10_SL1'
	# filOUT = path + 'out/IC10_SL1'

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
