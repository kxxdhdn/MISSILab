#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

PROCESS IMage

"""

import os
import math
import numpy as np
from scipy.io import readsav
from reproject import reproject_interp
import subprocess as SP

## astylo
from sinout import read_fits, write_fits, WCSextract, read_csv, write_csv, read_ascii
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


def wclean(filIN, cmod='eq', cfile=None, \
	wmod=0, filOUT=None, verbose=False):
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
	hdr, data, wave = read_fits(filIN)
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

		write_fits(filOUT, hdr, data_new, wave_new, wmod) # hdr auto changed
		
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
		'''
		self: filIN, wmod, hdr, w, is3d, Nx, Ny, Nw, im, wvl
		'''
		
		## INPUTS
		self.filIN = filIN
		self.wmod = wmod

		## read image/cube
		## self.hdr is a 2D (reduced) header
		self.hdr, self.w, self.is3d = WCSextract(filIN)
		self.Nx = self.hdr['NAXIS1']
		self.Ny = self.hdr['NAXIS2']
		print("Raw size (pix): {} * {}".format(
			self.Nx, self.Ny))
		## 3D cube slicing
		if self.is3d==True:
			self.im, self.wvl = read_fits(filIN)[1:3]
			self.Nw = len(self.wvl)
		else:
			self.im = read_fits(filIN)[1]
			self.wvl = None

	def addunc(self, uncIN=None):
		## uncertainty adding
		mu, sigma = 0., 1.

		if uncIN is not None:
			unc = read_fits(uncIN)[1]
			theta = np.random.normal(mu, sigma, unc.shape)
			## unc has the same dimension with im
			self.im += theta * unc

		return self.im

	def slice(self, filSL, suffix=None):
		## 3D cube slicing
		lstOUT = []
		if self.is3d==True:
			for k in range(self.Nw):
				## output filename list
				f = filSL+'_'+'0'*(4-len(str(k)))+str(k)
				if suffix is not None:
					f += suffix
				lstOUT.append(f)
				comment = "NO.{} image [SLICE]d from {}.fits".format(k, self.filIN)
				write_fits(f, self.hdr, self.im[k,:,:]) # addunc inclu
		else:
			print('Input file is a 2D image which cannot be sliced! ')
			f = filSL+'_0000'
			if suffix is not None:
				f += suffix
			lstOUT.append(f)
			write_fits(f, self.hdr, self.im) # addunc inclu
			print('Rewritten with only random noise added (if provided).')

		return lstOUT
	
	def crop(self, sizpix=None, cenpix=None, sizval=None, cenval=None, filOUT=None):
		'''
		If pix and val co-exist, pix will be taken.

		--- INPUT ---
		sizpix      crop size in pix (dx, dy)
		cenpix      crop center in pix (x, y)
		sizval      crop size in deg (dRA, dDEC) -> (dx, dy)
		cenval      crop center in deg (RA, DEC) -> (x, y)
		--- OUTPUT ---
		'''
		## Crop center
		if cenpix is None:
			## Convert coord
			cenpix = self.w.all_world2pix(cenval[0], cenval[1], 1)
		else:
			cenval = self.w.all_pix2world(np.array([cenpix]), 1)[0]
		if not (0<cenpix[0]<self.Nx and 0<cenpix[1]<self.Ny):
			print("ERROR: crop centre overpassed image border! ")
			exit()
		print('----------')
		print("Crop centre (RA, DEC): [{:.8}, {:.8}]".format(*cenval))
		print("Crop centre (x, y): [{}, {}]".format(*cenpix))
		
		## Crop size
		if sizpix is None:
			print("Crop size (dRA, dDEC): [{}, {}].".format(*sizval))
			## CDELTn needed (Physical increment at the reference pixel)
			sizpix = np.array(sizval) / np.array(self.hdr['CDELT1'], self.hdr['CDELT2'])
		print("Crop size (dx, dy): [{}, {}].".format(*sizpix))
		print('----------')
		
		## lowerleft origin
		xmin = math.floor(cenpix[0] - sizpix[0]/2.)
		ymin = math.floor(cenpix[1] - sizpix[1]/2.)
		xmax = xmin + sizpix[0]
		ymax = ymin + sizpix[1]
		if not (xmin>=0 and xmax<=self.Nx and ymin>=0 and ymax<=self.Ny):
			print("ERROR: crop region overpassed image border! ")
			exit()

		## OUTPUTS
		##---------
		## New image
		if self.is3d==True:
			self.im = self.im[:, ymin:ymax, xmin:xmax] # addunc inclu
			## recover 3D non-reduced header
			self.hdr = read_fits(self.filIN)[0]
		else:
			self.im = self.im[ymin:ymax, xmin:xmax] # addunc inclu
		## Modify header
		## Suppose no non-linear distortion
		self.hdr['CRPIX1'] = sizpix[0] / 2.
		self.hdr['CRPIX2'] = sizpix[1] / 2.
		self.hdr['CRVAL1'] = cenval[0]
		self.hdr['CRVAL2'] = cenval[1]
		## Write cropped image/cube
		if filOUT is not None:
			comment = "[ICROP]ped at centre: [{:.8}, {:.8}]. ".format(*cenval)
			# comment = "with size [{}, {}] (pix).".format(*sizpix)

			write_fits(filOUT, self.hdr, self.im, self.wvl, self.wmod, \
				COMMENT=comment)

		return self.im

class islice(improve):
	'''
	Slice a cube

	self: slist, path_tmp, (filIN, wmod, hdr, w, is3d, Nx, Ny, Nw, im, wvl)
	'''
	def __init__(self, filIN, filSL=None, uncIN=None, suffix=None):
		super().__init__(filIN)

		if filSL is None:
			path_tmp = os.getcwd()+'/tmp_slice/'
			if not os.path.exists(path_tmp):
				os.makedirs(path_tmp)
			self.path_tmp = path_tmp

			filSL = path_tmp+'slice'
		
		self.addunc(uncIN)

		self.slist = self.slice(filSL, suffix) # addunc inclu

	def image(self):
		return self.im

	def wave(self):
		return self.wvl

	def filenames(self):
		return self.slist

	def clean(self):
		fclean(path_tmp)

class icrop(improve):
	'''
	CROP 2D image or 3D cube
	'''
	def __init__(self, filIN, filOUT=None, \
		sizpix=None, cenpix=None, sizval=None, cenval=None, \
		uncIN=None, wmod=0):
		## slicrop: slice 
		super().__init__(filIN, wmod)

		self.addunc(uncIN)
		
		self.cropim = self.crop(sizpix=sizpix, cenpix=cenpix, \
			sizval=sizval, cenval=cenval, filOUT=filOUT) # addunc inclu

	def image(self):
		return self.cropim

	def wave(self):
		return self.wvl

class iproject(improve):
	'''
	PROJECT 2D image or 3D cube

	--- INPUT ---

	--- OUTPUT ---
	'''
	def __init__(self, filIN, filREF=None, hdREF=None, fmod='ref'):
		'''
		self: hdREF, (filIN, wmod, hdr, w, is3d, Nx, Ny, Nw, im, wvl)
		'''
		super().__init__(filIN)
		self.slist = None

		## Prepare reprojection header
		if filREF is not None:
			hdREF = WCSextract(filREF)[0]
			hdREF['EQUINOX'] = 2000.0

		if hdREF is not None:
			## Frame mode options
			##--------------------
			## Input WCS
			p = [[0, 0]]
			p.append([0, self.Ny])
			p.append([self.Nx, 0])
			p.append([self.Nx, self.Ny])
			val = self.w.all_pix2world(np.array(p), 1)
			## Ref WCS
			w = WCSextract(header=hdREF)[1]
			pts = []
			for i in range(4):
				pts.append(w.all_world2pix(val[i,0], val[i,1], 1))
			pts = np.array(pts)
			xmin = min(pts[:,0])
			xmax = max(pts[:,0])
			ymin = min(pts[:,1])
			ymax = max(pts[:,1])
			## Modify ref header
			if fmod=='ref':
				pass
			elif fmod=='rec': 
				hdREF['CRPIX1'] = -xmin
				hdREF['CRPIX2'] = -ymin
				hdREF['NAXIS1'] = math.ceil(xmax - xmin)
				hdREF['NAXIS2'] = math.ceil(ymax - ymin)
			elif fmod=='ext':
				if xmin<0:
					hdREF['CRPIX1'] = -xmin
				if ymin<0:
					hdREF['CRPIX2'] = -ymin
				hdREF['NAXIS1'] = math.ceil(max(xmax, hdREF['NAXIS1']-xmin, \
					xmax-xmin, hdREF['NAXIS1']))
				hdREF['NAXIS2'] = math.ceil(max(ymax, hdREF['NAXIS2']-ymin, \
					ymax-ymin, hdREF['NAXIS2']))
			## Save hdREF
			self.hdREF = hdREF
		else:
			print('ERROR: Can not find projection reference! ')
			exit()

	def reproject(self, filOUT=None, wmod=0, uncIN=None, filTMP=None):
		'''
		PROJECT 2D image or 3D cube

		--- INPUT ---
		fmod        output file frame
		            'ref' - same as ref frame (Default)
		            'rec' - recenter back to input frame
		            'ext' - cover both input and ref frame
		--- OUTPUT ---

		self: wmod, im, hdr, slist
		'''
		self.wmod = wmod
		self.addunc(uncIN)

		if filTMP is None:
			path_tmp = os.getcwd()+'/tmp_pro/'
			if not os.path.exists(path_tmp):
				os.makedirs(path_tmp)
			self.path_tmp = path_tmp

			filTMP = path_tmp+'slice'
		
		self.slist = self.slice(filTMP, '_') # addunc inclu

		## Do reprojection
		##-----------------
		im = []
		for f in self.slist:
			im.append(reproject_interp(f+'.fits', self.hdREF)[0])
			write_fits(f+'_pro', self.hdREF, im)
			fclean(f+'.fits')

		self.im = np.array(im)

		if filOUT is not None:
			self.hdr = self.hdREF

			comment = "Reprojected by <iproject>. "

			write_fits(filOUT, self.hdr, self.im, self.wvl, self.wmod, \
				COMMENT=comment)
	
		return self.im

	def wave(self):
		return self.wvl

	def filenames(self):
		return self.slist

	def clean(self):
		fclean(path_tmp)

class iconvolve(improve):
	'''
	(IDL based) CONVOLVE 2D image or 3D cube with given kernels

	--- INPUT ---
	filKER      convolution kernel(s) (tuple or list)
	--- OUTPUT ---
	'''
	def __init__(self, filIN, filKER, saveKER, \
		uncIN=None, psf=None, filTMP=None, wmod=0, filOUT=None):
		## INPUTS
		super().__init__(filIN, wmod)
		
		self.addunc(uncIN)

		## input kernel file list
		self.filKER = filKER
		## doc (csv) file of kernel list
		self.saveKER = saveKER
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
				f2conv = self.slice(self.filTMP) # addunc inclu
			else:
				f2conv = self.slice(self.filIN) # addunc inclu
			
			self.spitzer_irs()

		else:
			if self.filTMP is not None: # ?? TO DELETE
				f2conv = [self.filTMP]
			else:
				f2conv = [self.filIN]
			
		self.choker(f2conv)

		SP.call('cd '+ipath+'\nidl conv.pro', shell=True)

		## OUTPUTS
		##---------
		if self.is3d==True:
			im = []
			self.slist = []
			for f in f2conv:
				im.append(read_fits(f+'_conv')[1])
				self.slist.append(f+'_conv')

			self.convim = np.array(im)
			## recover non-reduced 3D header
			self.hdr = read_fits(self.filIN)[0]
		else:
			self.convim = read_fits(self.filIN+'_conv')[1]
		
		if self.filOUT is not None:
			comment = "Convolved by G. Aniano's IDL routine."
			# comment = 'https://www.astro.princeton.edu/~ganiano/Kernels.html'
			write_fits(self.filOUT, self.hdr, self.convim, self.wvl, self.wmod, \
				COMMENT=comment)

	def image(self):
		return self.convim

	def wave(self):
		return self.wvl

	def filenames(self):
		return self.slist

class slitextract(improve):
	'''
	AKARI/IRC spectroscopy slit coord extraction

	--- INPUT ---
	filOUT      output FITS file
	dirIRC      path of IRC dataset
	parobs[0]   observation id
	parobs[1]   IRC N3 (long exp) frame (2MASS corrected; 90 deg rot needed)
	parobs[2]   slit
	Nw          num of wave
	Ny          slit length
	Nx          slit width
	--- OUTPUT ---
	'''
	def __init__(self, dirIRC, parobs, Nw=259, Ny=31, Nx=None, filOUT=None):
		path = dirIRC + parobs[0] + '/irc_specred_out_' + parobs[2]+'/'
		filSAV = path + parobs[0] + '.N3_NG.IRC_SPECRED_OUT'
		filIN = path + parobs[1]
		super().__init__(filIN)
		
		## Read output data
		spec_arr = []
		for i in range(Ny):
			spec_arr.append(read_ascii(path+'spec'+str(i), float, '.spc'))
		spec_arr = np.array(spec_arr)
		table = readsav(filSAV+'.sav', python_dict=True)['source_table']
		ref_x = table['image_y'][0] # slit ref x
		ref_y = 512-table['image_x'][0] # slit ref y
		## Slit width will be corrected by intercalib with IRS data after reprojection
		if parobs[2]=='Ns':
			Nx = 3 # 412pix * 5"/10' = 3.43 pix (Ns)
		elif parobs[2]=='Nh':
			Nx = 2 # 412pix * 3"/10' = 2.06 pix (Nh)
		
		cube = np.empty([Nw,Ny,Nx])
		unc = np.empty([Nw,Ny,Nx])
		wave = np.empty(Nw)
		
		## Alternative (see IRC_SPEC_TOOL; not needed if use IDL extracted spec)
		'''
		image = readsav(savIN+'.sav', python_dict=True)['specimage_n_wc']
		image = image[::-1] # -> ascending order
		noise = readsav(savIN+'.sav', python_dict=True)['noisemap_n']
		noise = noise[::-1]
		wave = readsav(savIN+'.sav', python_dict=True)['wave_array']
		wave = wave[::-1] # -> ascending order
		Nw = image.shape[0] # num of wave
		Ny = image.shape[1] # slit length
		spec_y = table['spec_y'][0] # ref pts of wavelength
		
		d_wave_offset_pix = -(spec_y-round(spec_y[0])) # Wave shift
		warr = np.arange(Nw)
		wave_shift = np.interp(warr+d_wave_offset_pix, warr, wave)
		
		for k in range(Nw):
			for j in range(Ny):
				for i in range(Nx):
					cube[k][j][i] = image[k][j]
					unc[k][j][i] = noise[k][j]
		'''
		
		## Broaden cube width
		for k in range(Nw):
			for j in range(Ny):
				for i in range(Nx):
					cube[k][j][i] = spec_arr[j,k,1]
					unc[k][j][i] = (spec_arr[j,k,3]-spec_arr[j,k,2])/2
			wave[k] = spec_arr[0,k,0]

		## Save spec in wave ascending order
		self.slitim = cube[::-1]
		self.unc = unc[::-1]
		self.wvl = wave[::-1]
		
		## Modify cube header
		self.crop(sizpix=(Nx, Ny), cenpix=(ref_x, ref_y))
		# self.hdr['CTYPE3'] = 'WAVE-TAB'
		self.hdr['CUNIT1'] = 'deg'
		self.hdr['CUNIT2'] = 'deg'
		self.hdr['BUNIT'] = 'MJy/sr'
		self.hdr['EQUINOX'] = 2000.0
		comment = "Assembled AKARI/IRC slit spectroscopy cube. "
		uncom = "Assembled AKARI/IRC slit spectroscopy uncertainty cube. "

		write_fits(filOUT, self.hdr, self.slitim, self.wvl, \
			COMMENT=comment)
		write_fits(filOUT+'_unc', self.hdr, self.unc, self.wvl, \
			COMMENT=uncom)
		print('obs_id: {} \nslit: {}'.format(parobs[0], parobs[2]))
		print('----------\n')

	def image(self):
		return self.slitim

	def wave(self):
		return self.wvl

class imontage(iproject):
	'''
	<improve> based montage tool

	--- INPUT ---
	flist       list of files to combine
	--- OUTPUT ---
	'''
	def __init__(self, flist, ref=None, hdr_ref=None, fmod='ext'):

		super().__init__(flist[0], ref, hdr_ref, fmod)
		self.flist = flist
		self.fmod = fmod
	
	def repro(self, out=None, wmod=0, unc=None, tmp=None):
		'''
		Reproject input file
		'''
		image = self.reproject(out, wmod, unc, tmp)

		return image

	# def pointing(self, image):

	def combine(self, method='average', template=None, filOUT=None, ulist=None):
		'''
		Combine the reprojected file to the new input file (after reproj it)

		'''
		Nf = len(self.flist)
		template = os.getcwd()+'/repro_template'

		## Extend reprojection header
		for i in range(Nf):
			## self.hdREF renewed in every cycle
			super().__init__(self.flist[i], None, self.hdREF, self.fmod)
			if i<Nf-1:
				self.repro() # unc will be added later
			else:
				self.repro(out=template) # save template

		## Reproject images
		hyperim = []
		for i in range(Nf):
			super().__init__(self.flist[i], None, self.hdREF) # fmod='ref'
			hyperim.append(self.repro(unc=ulist[i]))
		hyperim = np.array(hyperim)

		# ## Alternative reprojection (with repro_template.fits)
		# hyperim = []
		# for i in range(Nf):
		# 	super().__init__(self.flist[i], template) # fmod='ref'
		# 	hyperim.append(self.repro(unc=ulist[i]))
		# hyperim = np.array(hyperim)

		## Combine images
		Nw, Ny, Nx = hyperim[0].shape
		combim = np.zeros((Nw, Ny, Nx))
		for k in range(Nw):
			for y in range(Ny):
				for x in range(Nx):
					combim[k,y,x] = np.nanmean(hyperim[:,k,y,x])

		if filOUT is not None:
			comment = "A <imontage> production"

			write_fits(filOUT, self.hdREF, combim, self.wvl, \
				COMMENT=comment)

		# self.clean()
		
"""
------------------------------ MAIN (test) ------------------------------
"""
if __name__ == "__main__":

	from myfunclib import MCerror

	import os
	path_test = os.getcwd() + '/../tests/'

	## Test slitextract + imontage
	##-----------------------------
	path_data = '/Users/dhu/Data/AKARI/data/'

	obs_id = ['1420415.1', '1420415.2']
	N3 = ['F011213824_N002', 'F011213865_N002']
	slit = ['Ns', 'Nh']

	par_obs = []
	file_out = []
	for i, obs in enumerate(obs_id):
		for s in slit:
			par_obs.append([obs, N3[i], s])
			file_out.append(path_test + 'out/M83_' + obs + '_' + s)

	spec = []
	unc_out = []
	for j, par in enumerate(par_obs):
		spec.append(slitextract(path_data, par, filOUT=file_out[j]))
		unc_out.append(file_out[j]+'_unc')
	# print(spec[0].image().shape, spec[0].wave())
	
	mont = imontage(file_out, file_out[0])
	mont.combine(filOUT=path_test + 'out/M83_IRC', ulist=unc_out)

	## Test FITS ref point shift (improve.crop)
	##-----------------------------------------------
	# filIN = path_test + 'data/M83_0'
	# filOUT = path_test + 'out/M83_ref_shift'
	
	# w = WCSextract(filIN)[1]
	# hdr, data, wave = read_fits(filIN)
	# hdr['CRPIX1'] = hdr['NAXIS1'] / 2.
	# hdr['CRPIX2'] = hdr['NAXIS2'] / 2.
	# pix = np.array([(hdr['CRPIX1'], hdr['CRPIX2'])])
	# hdr['CRVAL1'], hdr['CRVAL2'] = w.all_pix2world(pix, 1)[0]

	# write_fits(filOUT, hdr, data, wave)

	## Test interfill
	##----------------
	# Nmc = 2

	# filIN = path_test + 'data/IC10_SL1'
	# filOUT = path_test + 'out/IC10_SL1'

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
