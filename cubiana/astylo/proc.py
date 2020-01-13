#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

image PROCessing

"""

import os
import math
import numpy as np
from scipy.io import readsav
from reproject import reproject_interp
import subprocess as SP

## Local
from bio import read_fits, write_fits, ext_wcs, read_csv, write_csv, read_ascii
from lib import fclean, closest, bsplinterpol

savext = '.sav'

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
	ds = read_fits(filIN)
	hdr = ds.header
	data = ds.data
	wave = ds.wave
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
	ds = read_fits(filIN)
	oldim = ds.data
	hdr = ds.header
	w = ext_wcs(filIN).WCS
	# hdr['NAXIS1'] = x1 - x0 + 1
	# hdr['NAXIS2'] = y1 - y0 + 1
	hdr['CRPIX1'] += -x0
	hdr['CRPIX2'] += -y0
	newim = oldim[y0:y1+1, x0:x1+1]

	write_fits(filOUT, hdr, newim)

	return newim

##-----------------------------------------------

##			<impro> based tools

##-----------------------------------------------

class impro:
	'''
	IMage PROcessing VEssel
	'''
	def __init__(self, filIN, wmod=0, verbose=True):
		'''
		self: filIN, wmod, hdr, w, dim, Nx, Ny, Nw, im, wvl
		'''
		
		## INPUTS
		self.filIN = filIN
		self.wmod = wmod
		self.verbose = verbose

		## read image/cube
		## self.hdr is a 2D (reduced) header
		ws = ext_wcs(filIN)
		self.hdr = ws.header
		self.w = ws.WCS
		self.Nx = self.hdr['NAXIS1']
		self.Ny = self.hdr['NAXIS2']
		self.Nw = None
		if verbose==True:
			print('<impro> file: ', filIN)
			print('Raw size (pix): {} * {}'.format(self.Nx, self.Ny))
		## 3D cube slicing
		ds = read_fits(filIN)
		self.im = ds.data
		self.wvl = ds.wave
		self.dim = len(self.im.shape)
		if self.dim==3:
			self.Nw = len(self.wvl)

	def rand_norm(self, file=None, unc=None, sigma=1., mu=0.):
		'''
		Add random N(0,1) noise
		'''
		if file is not None:
			unc = read_fits(file).data

		if unc is not None:
			## unc should have the same dimension with im
			ax = unc.shape
			theta = np.random.normal(mu, sigma, ax)
			self.im += theta * unc

		return self.im

	def rand_splitnorm(self, file=None, unc=None, sigma=1., mu=0.):
		'''
		Add random SN(0,lam,lam*tau) noise

		------ INPUT ------
		file                2 FITS files for unc of left & right sides
		unc                 2 uncertainty ndarrays
		------ OUTPUT ------
		'''
		if file is not None:
			unc = []
			for f in file:
				unc.append(read_fits(f).data)
			
		if unc is not None:
			## unc[i] should have the same dimension with self.im
			tau = unc[1]/unc[0]
			peak = 1/(1+tau)
			ax = tau.shape
			theta = np.random.normal(mu, sigma, ax) # ~N(0,1)
			flag = np.random.random(ax) # ~U(0,1)
			if len(ax)==2:
				for x in range(ax[1]):
					for y in range(ax[0]):
						if flag[y,x]<peak[y,x]:
							self.im[y,x] += -abs(theta[y,x]) * unc[0][y,x]
						else:
							self.im[y,x] += abs(theta[y,x]) * unc[1][y,x]
			elif len(ax)==3:
				for x in range(ax[2]):
					for y in range(ax[1]):
						for k in range(ax[0]):
							if flag[k,y,x]<peak[k,y,x]:
								self.im[k,y,x] += -abs(theta[k,y,x]) * unc[0][k,y,x]
							else:
								self.im[k,y,x] += abs(theta[k,y,x]) * unc[1][k,y,x]

		return self.im

	def slice(self, filSL, suffix=None):
		## 3D cube slicing
		slist = []
		if self.dim==3:
			for k in range(self.Nw):
				## output filename list
				f = filSL+'_'+'0'*(4-len(str(k)))+str(k)
				if suffix is not None:
					f += suffix
				slist.append(f)
				comment = "NO.{} image [SLICE]d from {}.fits".format(k, self.filIN)
				write_fits(f, self.hdr, self.im[k,:,:]) # gauss_noise inclu
		else:
			print('Input file is a 2D image which cannot be sliced! ')
			f = filSL+'_0000'
			if suffix is not None:
				f += suffix
			slist.append(f)
			write_fits(f, self.hdr, self.im) # gauss_noise inclu
			print('Rewritten with only random noise added (if provided).')

		return slist
	
	def crop(self, filOUT=None, sizpix=None, cenpix=None, sizval=None, cenval=None):
		'''
		If pix and val co-exist, pix will be taken.

		------ INPUT ------
		filOUT              output file
		sizpix              crop size in pix (dx, dy)
		cenpix              crop center in pix (x, y)
		sizval              crop size in deg (dRA, dDEC) -> (dx, dy)
		cenval              crop center in deg (RA, DEC) -> (x, y)
		------ OUTPUT ------
		self.im             cropped image array
		'''
		## Crop center
		##-------------
		if cenpix is None:
			if cenval is None:
				print('ERROR: crop center unavailable! ')
				exit()
			else:
				## Convert coord
				cenpix = self.w.all_world2pix(cenval[0], cenval[1], 1)
		else:
			cenval = self.w.all_pix2world(np.array([cenpix]), 1)[0]
		if not (0<cenpix[0]<self.Nx and 0<cenpix[1]<self.Ny):
			print("ERROR: crop centre overpassed image border! ")
			exit()

		## Crop size
		##-----------
		if sizpix is None:
			if sizval is None:
				print('ERROR: crop size unavailable! ')
				exit()
			else:
				## CDELTn needed (Physical increment at the reference pixel)
				sizpix = np.array(sizval) / np.array(self.hdr['CDELT1'], self.hdr['CDELT2'])
		# else:
			# sizval = np.array(sizpix) * np.array(self.hdr['CDELT1'], self.hdr['CDELT2'])

		if self.verbose==True:
			print('----------')
			print("Crop centre (RA, DEC): [{:.8}, {:.8}]".format(*cenval))
			# print("Crop size (dRA, dDEC): [{}, {}]\n".format(*sizval))
			print("Crop centre (x, y): [{}, {}]".format(*cenpix))
			print("Crop size (dx, dy): [{}, {}]".format(*sizpix))
			print('----------')
		
		## Lowerleft origin
		##------------------
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
		if self.dim==3:
			self.im = self.im[:, ymin:ymax, xmin:xmax] # gauss_noise inclu
			## recover 3D non-reduced header
			self.hdr = read_fits(self.filIN).header
		elif self.dim==2:
			self.im = self.im[ymin:ymax, xmin:xmax] # gauss_noise inclu
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

class islice(impro):
	'''
	Slice a cube

	self: slist, path_tmp, (filIN, wmod, hdr, w, dim, Nx, Ny, Nw, im, wvl)
	'''
	def __init__(self, filIN, filSL=None, filUNC=None, dist='norm', suffix=None):
		super().__init__(filIN)

		if filSL is None:
			cur_tmp = 1 # slices stocked in current path
			path_tmp = os.getcwd()+'/tmp_proc/'
			if not os.path.exists(path_tmp):
				os.makedirs(path_tmp)
			self.path_tmp = path_tmp

			filSL = path_tmp+'slice'
		
		if dist=='norm':
			self.rand_norm(filUNC)
		elif dist=='splitnorm':
			self.rand_splitnorm(filUNC)

		self.slist = self.slice(filSL, suffix) # gauss_noise inclu

	def image(self):
		return self.im

	def wave(self):
		return self.wvl

	def filenames(self):
		return self.slist

	def clean(self, file=None):
		if file is not None:
			fclean(file)
		else:
			fclean(self.path_tmp)

class icrop(impro):
	'''
	CROP 2D image or 3D cube
	'''
	def __init__(self, filIN, filOUT=None, \
		sizpix=None, cenpix=None, sizval=None, cenval=None, \
		filUNC=None, dist='norm', wmod=0):
		## slicrop: slice 
		super().__init__(filIN, wmod)

		if dist=='norm':
			self.rand_norm(filUNC)
		elif dist=='splitnorm':
			self.rand_splitnorm(filUNC)
		
		im_crop = self.crop(filOUT=filOUT, sizpix=sizpix, cenpix=cenpix, \
			sizval=sizval, cenval=cenval) # gauss_noise inclu

	def image(self):
		return self.im

	def wave(self):
		return self.wvl

class iproject(impro):
	'''
	PROJECT 2D image or 3D cube
	i means <impro>-based or initialize

	------ INPUT ------
	filIN               input FITS file
	filREF              ref file (priority if co-exist with input header)
	hdREF               ref header
	fmod                output file frame
		                  'ref' - same as ref frame (Default)
		                  'rec' - recenter back to input frame
		                  'ext' - cover both input and ref frame
	------ OUTPUT ------
	'''
	def __init__(self, filIN, filREF=None, hdREF=None, fmod='ref'):
		'''
		self: hdr_ref, (filIN, wmod, hdr, w, dim, Nx, Ny, Nw, im, wvl)
		'''
		super().__init__(filIN)

		## Set path of tmp files
		path_tmp = os.getcwd()+'/tmp_proc/'
		if not os.path.exists(path_tmp):
			os.makedirs(path_tmp)
		self.path_tmp = path_tmp

		## Prepare reprojection header
		if filREF is not None:
			hdREF = ext_wcs(filREF).header
			hdREF['EQUINOX'] = 2000.0

		if hdREF is not None:
			## Frame mode (fmod) options
			##---------------------------
			if fmod=='ref':
				pass
			else:
				## Input WCS (old)
				pix_old = [[0, 0]]
				pix_old.append([0, self.Ny])
				pix_old.append([self.Nx, 0])
				pix_old.append([self.Nx, self.Ny])
				world_arr = self.w.all_pix2world(np.array(pix_old), 1)
				## Ref WCS (new)
				w = ext_wcs(header=hdREF).WCS
				pix_new = w.all_world2pix(world_arr, 1)
				xmin = min(pix_new[:][0])
				xmax = max(pix_new[:,0])
				ymin = min(pix_new[:,1])
				ymax = max(pix_new[:,1])

				## Modify ref header
				if fmod=='rec': 
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
			self.hdr_ref = hdREF
		else:
			print('ERROR: Can not find projection reference! ')
			exit()

	def reproject(self, filOUT=None, filUNC=None, dist='norm', wmod=0, filTMP=None):
		'''
		PROJECT 2D image or 3D cube

		------ INPUT ------
		filOUT              output FITS file
		wmod                wave mode of output file
		filUNC              unc file (add gaussian noise)
		filTMP              save tmp files
		------ OUTPUT ------

		self: filTMP
		'''
		self.wmod = wmod
		if dist=='norm':
			self.rand_norm(filUNC)
		elif dist=='splitnorm':
			self.rand_splitnorm(filUNC)

		if filTMP is None:
			cur_tmp = 1 # tmp files stocked in current path
			filTMP = self.path_tmp+'slice'
		else:
			self.filTMP = filTMP
		
		self.slist = self.slice(filTMP, '_') # gauss_noise inclu

		## Do reprojection
		##-----------------
		cube_rep = []
		for f in self.slist:
			im_rep = reproject_interp(f+'.fits', self.hdr_ref)[0]
			cube_rep.append(im_rep)

			write_fits(f+'_rep', self.hdr_ref, im_rep)
			fclean(f+'.fits')
		if cur_tmp==1:
			fclean(self.path_tmp)

		self.im = np.array(cube_rep)

		if filOUT is not None:
			self.hdr = self.hdr_ref

			comment = "Reprojected by <iproject>. "
			write_fits(filOUT, self.hdr, self.im, self.wvl, self.wmod, \
				COMMENT=comment)
	
		return self.im

	def wave(self):
		return self.wvl

	def filenames(self):
		return self.slist

	def clean(self, file=None):
		if file is not None:
			fclean(file)
		else:
			fclean(self.path_tmp)

class imontage(iproject):
	'''
	<impro> based montage tool

	------ INPUT ------
	flist               list of files
	ref                 
	header_ref          
	fmod                
	------ OUTPUT ------
	'''
	def __init__(self, flist, ref=None, header_ref=None, fmod='ext'):
		'''
		self: hdr_ref, path_tmp, (filIN, wmod, hdr, w, dim, Nx, Ny, Nw, im, wvl)
		'''
		super().__init__(flist[0], ref, header_ref, fmod)
		self.flist = flist

		## Set path of tmp files
		path_tmp = os.getcwd()+'/tmp_proc/'
		if not os.path.exists(path_tmp):
			os.makedirs(path_tmp)
		self.path_tmp = path_tmp

		## Refresh self.hdr_ref in every circle
		if fmod=='ext':
			for f in flist:
				super().__init__(f, None, self.hdr_ref, 'ext')

	def footprint(self, filOUT=None, wmod=0):
		'''
		self: hdr_ref, (filIN, wmod, hdr, w, dim, Nx, Ny, Nw, im, wvl)
		'''
		if filOUT is None:
			cur_tmp = 1 # tmp files stocked in current path
			filOUT = self.path_tmp+'footprint'
		
		Nx = self.hdr_ref['NAXIS1']
		Ny = self.hdr_ref['NAXIS2']
		im_fp = np.ones((Ny, Nx))
		
		comment = "<imontage> footprint"
		write_fits(filOUT, self.hdr_ref, im_fp, None, wmod, \
				COMMENT=comment)

		return im_fp

	def combine(self, filOUT=None, method='average', ulist=None, dist='norm'):
		'''
		Reproject and combine input files to the ref WCS

		'''
		## Reproject images
		hyperim = []
		for i, f in enumerate(self.flist):
			super().__init__(f, None, self.hdr_ref) # fmod='ref'
			hyperim.append(self.reproject(filUNC=ulist[i], dist=dist))
		hyperim = np.array(hyperim)

		## Combine images
		# dim = len(hyperim[0].shape)
		if self.dim==2:
			Ny, Nx = hyperim[0].shape
			im_comb = np.zeros((Ny, Nx))
			for y in range(Ny):
				for x in range(Nx):
					im_comb[y,x] = np.nanmean(hyperim[:,y,x])
		elif self.dim==3:
			Nw, Ny, Nx = hyperim[0].shape
			print(hyperim.shape)
			im_comb = np.zeros((Nw, Ny, Nx))
			for k in range(Nw):
				for y in range(Ny):
					for x in range(Nx):
						im_comb[k,y,x] = np.nanmean(hyperim[:,k,y,x])

		if filOUT is not None:
			comment = "A <imontage> production"
			write_fits(filOUT, self.hdr_ref, im_comb, self.wvl, \
				COMMENT=comment)

		return im_comb

	def clean(self, file=None):
		if file is not None:
			fclean(file)
		else:
			fclean(self.path_tmp)

class iconvolve(impro):
	'''
	Convolve 2D image or 3D cube with given kernels
	i means <impro>-based or IDL-based

	------ INPUT ------
	filIN               input FITS file
	filKER              convolution kernel(s) (tuple or list)
	kfile               CSV file stocking kernel names
	filUNC              unc file (add gaussian noise)
	psf                 PSF list
	filTMP              temporary file
	wmod                wave mode
	filOUT              output file
	------ OUTPUT ------
	'''
	def __init__(self, filIN, filKER, kfile, \
		filUNC=None, dist='norm', psf=None, filTMP=None, wmod=0, filOUT=None):
		## INPUTS
		super().__init__(filIN, wmod)
		
		if dist=='norm':
			self.rand_norm(filUNC)
		elif dist=='splitnorm':
			self.rand_splitnorm(filUNC)

		## input kernel file list
		self.filKER = filKER
		## doc (csv) file of kernel list
		self.kfile = kfile
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
		write_csv(self.kfile, header=['Images', 'Kernels'], dataset=klist)

	def do_conv(self, ipath):
		'''
		------ INPUT ------
		ipath               path of IDL routines
		------ OUTPUT ------
		'''
		if self.dim==3:
			if self.filTMP is not None:
				f2conv = self.slice(self.filTMP) # gauss_noise inclu
			else:
				f2conv = self.slice(self.filIN) # gauss_noise inclu
			
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
		if self.dim==3:
			im = []
			self.slist = []
			for f in f2conv:
				im.append(read_fits(f+'_conv').data)
				self.slist.append(f+'_conv')

			self.convim = np.array(im)
			## recover non-reduced 3D header
			self.hdr = read_fits(self.filIN).header
		else:
			self.convim = read_fits(self.filIN+'_conv').data
		
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

class sextract(impro):
	'''
	AKARI/IRC spectroscopy slit coord extraction
	s means slit, spectral cube or SAV file

	------ INPUT ------
	filOUT              output FITS file
	ipath               path of IRC dataset
	parobs[0]           observation id
	parobs[1]           slit name
	parobs[2]           IRC N3 (long exp) frame (2MASS corrected; 90 deg rot needed)
	Nw                  num of wave
	Ny                  slit length
	Nx                  slit width
	------ OUTPUT ------
	'''
	def __init__(self, ipath=None, parobs=None):
		self.path = ipath + parobs[0] + '/irc_specred_out_' + parobs[1]+'/'
		filIN = self.path + parobs[2]
		super().__init__(filIN)
		
		self.filSAV = self.path + parobs[0] + '.N3_NG.IRC_SPECRED_OUT'
		self.table = readsav(self.filSAV+savext, python_dict=True)['source_table']

		## Slit width will be corrected by intercalib with IRS data after reprojection
		if parobs[1]=='Ns':
			self.slit_width = 3 # 412pix * 5"/10' = 3.43 pix (Ns)
		elif parobs[1]=='Nh':
			self.slit_width = 2 # 412pix * 3"/10' = 2.06 pix (Nh)

		print('\n----------')
		print('Slit extracted from ')
		print('obs_id: {} \nslit: {}'.format(parobs[0], parobs[1]))
		print('----------\n')

	def rand_pointing(self, sigma=0.):
		'''
		Add pointing uncertainty to WCS

		------ INPUT ------
		sigma               pointing accuracy (deg)
		------ OUTPUT ------
		'''
		d_ro = abs(np.random.normal(0., sigma)) # N(0,sigma)
		d_phi = np.random.random() *2. * np.pi # U(0,2pi)
		self.hdr['CRVAL1'] += d_ro * np.cos(d_phi)
		self.hdr['CRVAL2'] += d_ro * np.sin(d_phi)

		return d_ro, d_phi

	def spec_build(self, filOUT=None, Ny=31, sig_pt=0.):
		'''
		Build the spectral cube/slit from spectra extracted by IDL pipeline
		(see IRC_SPEC_TOOL, plot_spec_with_image)
		'''
		Nx = self.slit_width
		ref_x = self.table['image_y'][0] # slit ref x
		ref_y = 512 - self.table['image_x'][0] # slit ref y

		## Get slit coord from 2MASS corrected N3 frame
		## Do NOT touch self.im (N3 frame, 2D) before this step
		self.crop(sizpix=(Nx, Ny), cenpix=(ref_x, ref_y))
		# self.hdr['CTYPE3'] = 'WAVE-TAB'
		self.hdr['CUNIT1'] = 'deg'
		self.hdr['CUNIT2'] = 'deg'
		self.hdr['BUNIT'] = 'MJy/sr'
		self.hdr['EQUINOX'] = 2000.0

		## Add pointing unc
		self.rand_pointing(sig_pt)

		## Read spec
		spec_arr = []
		for i in range(Ny):
			spec_arr.append(read_ascii(self.path+'spec'+str(i), float, '.spc'))
		spec_arr = np.array(spec_arr)
		Nw = len(spec_arr[0,:,0])
		
		## Broaden cube width
		cube = np.empty([Nw,Ny,Nx])
		unc = np.empty([Nw,Ny,Nx]) # Symmetric unc
		unc_N = np.empty([Nw,Ny,Nx]) # Asymmetric negtive
		unc_P = np.empty([Nw,Ny,Nx]) # Asymmetric positive
		wave = np.empty(Nw)
		for k in range(Nw):
			for j in range(Ny):
				for i in range(Nx):
					cube[k][j][i] = spec_arr[j,k,1]
					unc[k][j][i] = (spec_arr[j,k,3]-spec_arr[j,k,2])/2
					unc_N[k][j][i] = (spec_arr[j,k,1]-spec_arr[j,k,2])
					unc_P[k][j][i] = (spec_arr[j,k,3]-spec_arr[j,k,1])
			wave[k] = spec_arr[0,k,0]

		## Save spec in wave ascending order
		self.cube = cube[::-1]
		self.unc = unc[::-1]
		self.unc_N = unc_N[::-1]
		self.unc_P = unc_P[::-1]
		self.wvl = wave[::-1]

		if filOUT is not None:
			comment = "Assembled AKARI/IRC slit spectroscopy cube. "
			write_fits(filOUT, self.hdr, self.cube, self.wvl, \
				COMMENT=comment)

			uncom = "Assembled AKARI/IRC slit spec uncertainty cube. "
			write_fits(filOUT+'_unc', self.hdr, self.unc, self.wvl, \
				COMMENT=uncom)

			uncom_N = "Assembled AKARI/IRC slit spec uncertainty (N) cube. "
			write_fits(filOUT+'_unc_N', self.hdr, self.unc_N, self.wvl, \
				COMMENT=uncom)

			uncom_P = "Assembled AKARI/IRC slit spec uncertainty (P) cube. "
			write_fits(filOUT+'_unc_P', self.hdr, self.unc_P, self.wvl, \
				COMMENT=uncom)

		return self.cube

	def sav_build(self):
		'''
		Alternative extraction from SAV file
		Including wave calib, ?, etc. 
		(see IRC_SPEC_TOOL, plot_spec_with_image)
		'''
		filSAV = self.filSAV
		table = self.table
		## Read SAV file
		image = readsav(filSAV+savext, python_dict=True)['specimage_n_wc']
		image = image[::-1] # -> ascending order
		noise = readsav(filSAV+savext, python_dict=True)['noisemap_n']
		noise = noise[::-1]
		wave = readsav(filSAV+savext, python_dict=True)['wave_array']
		wave = wave[::-1] # -> ascending order
		Nw = image.shape[0] # num of wave
		Ny = image.shape[1] # slit length
		ref_x = table['image_y'][0] # slit ref x
		ref_y = 512-table['image_x'][0] # slit ref y
		spec_y = table['spec_y'][0] # ref pts of wavelength
		
		d_wave_offset_pix = -(spec_y-round(spec_y[0])) # Wave shift
		warr = np.arange(Nw)
		wave_shift = np.interp(warr+d_wave_offset_pix, warr, wave)
		
		for k in range(Nw):
			for j in range(Ny):
				for i in range(Nx):
					cube[k][j][i] = image[k][j]
					unc[k][j][i] = noise[k][j]

	def image(self):
		return self.cube

	def wave(self):
		return self.wvl
	
"""
------------------------------ MAIN (test) ------------------------------
"""
if __name__ == "__main__":

	import os
	path_test = os.getcwd() + '/../tests/'
