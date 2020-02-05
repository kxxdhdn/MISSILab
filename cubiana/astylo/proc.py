#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

image PROCessing

"""
from tqdm import tqdm, trange
import os
import math
import numpy as np
from scipy.io import readsav
from astropy import wcs
from reproject import reproject_interp
import subprocess as SP

## Local
from bio import read_fits, write_fits, read_csv, write_csv, read_ascii
from astrolib import fixwcs, pc2cd
from lib import fclean, nanavg, closest, bsplinterpol

savext = '.sav'


##-----------------------------------------------

##			<impro> based tools

##-----------------------------------------------

class impro:
	'''
	IMage PROcessing VEssel
	'''
	def __init__(self, filIN, wmod=0, verbose=False):
		'''
		self: filIN, wmod, hdr, w, dim, Nx, Ny, Nw, im, wvl
		'''
		
		## INPUTS
		self.filIN = filIN
		self.wmod = wmod
		self.verbose = verbose

		## read image/cube
		## self.hdr is a 2D (reduced) header
		ws = fixwcs(filIN)
		self.hdr = ws.header
		self.w = ws.wcs
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
			theta = np.random.normal(mu, sigma, self.im.shape)
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
			theta = np.random.normal(mu, sigma, self.im.shape) # ~N(0,1)
			flag = np.random.random(self.im.shape) # ~U(0,1)
			if self.dim==2:
				for x in range(self.Nx):
					for y in range(self.Ny):
						if flag[y,x]<peak[y,x]:
							self.im[y,x] += -abs(theta[y,x]) * unc[0][y,x]
						else:
							self.im[y,x] += abs(theta[y,x]) * unc[1][y,x]
			elif self.dim==3:
				for x in range(self.Nx):
					for y in range(self.Ny):
						for k in range(self.Nw):
							if flag[k,y,x]<peak[k,y,x]:
								self.im[k,y,x] += -abs(
									theta[k,y,x]) * unc[0][k,y,x]
							else:
								self.im[k,y,x] += abs(
									theta[k,y,x]) * unc[1][k,y,x]

		return self.im

	def slice(self, filSL, postfix=''):
		## 3D cube slicing
		slist = []
		if self.dim==3:
			hdr = self.hdr.copy()
			for kw in self.hdr.keys():
				if '3' in kw:
					del hdr[kw]
			hdr['NAXIS'] = 2
			for k in range(self.Nw):
				## output filename list
				f = filSL+'_'+'0'*(4-len(str(k)))+str(k)+postfix
				slist.append(f)
				write_fits(f, hdr, self.im[k,:,:]) # gauss_noise inclu
		else:
			print('Input file is a 2D image which cannot be sliced! ')
			f = filSL+'_0000'+postfix
			slist.append(f)
			write_fits(f, self.hdr, self.im) # gauss_noise inclu
			print('Rewritten with only random noise added (if provided).')

		return slist

	def slice_inv_sqrt(self, filSL, postfix=''):
		## Inversed square cube slicing
		inv_sqrt = 1./self.im**2
		slist = []
		if self.dim==3:
			hdr = self.hdr.copy()
			for kw in self.hdr.keys():
				if '3' in kw:
					del hdr[kw]
			hdr['NAXIS'] = 2
			for k in range(self.Nw):
				## output filename list
				f = filSL+'_'+'0'*(4-len(str(k)))+str(k)+postfix
				slist.append(f)
				write_fits(f, hdr, inv_sqrt[k,:,:]) # gauss_noise inclu
		else:
			f = filSL+'_0000'+postfix
			slist.append(f)
			write_fits(f, self.hdr, inv_sqrt) # gauss_noise inclu
		
		print('<impro> Inversed square cube slicing.')

		return slist
	
	def crop(self, filOUT=None, \
		sizpix=None, cenpix=None, sizval=None, cenval=None):
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
				try:
					cenpix = self.w.all_world2pix(cenval[0], cenval[1], 1)
				except wcs.wcs.NoConvergence as e:
					cenpix = e.best_solution
					print("Best solution:\n{0}".format(e.best_solution))
					print("Achieved accuracy:\n{0}".format(e.accuracy))
					print("Number of iterations:\n{0}".format(e.niter))
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
				sizpix = np.array(sizval) / np.array(
					self.hdr['CDELT1'], self.hdr['CDELT2'])
		# else:
			# sizval = np.array(sizpix) * np.array(
			# 	self.hdr['CDELT1'], self.hdr['CDELT2'])

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
	def __init__(self, filIN, filSL=None, filUNC=None, dist='norm', \
		slicetype=None, postfix=''):
		super().__init__(filIN)

		if filSL is None:
			path_tmp = os.getcwd()+'/tmp_proc/'
			if not os.path.exists(path_tmp):
				os.makedirs(path_tmp)
			self.path_tmp = path_tmp

			filSL = path_tmp+'slice'
		
		if dist=='norm':
			self.rand_norm(filUNC)
		elif dist=='splitnorm':
			self.rand_splitnorm(filUNC)

		if slicetype is None:
			self.slist = self.slice(filSL, postfix) # gauss_noise inclu
		elif slicetype=='inv_sqrt':
			self.slist = self.slice_inv_sqrt(filSL, postfix)

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

class imontage(impro):
	'''
	2D image or 3D cube montage toolkit
	i means <impro>-based or initialize

	------ INPUT ------
	file                FITS file (list, cf impro.filIN)
	filREF              ref file (priority if co-exist with input header)
	hdREF               ref header
	fmod                output image frame mode
		                  'ref' - same as ref frame (Default)
		                  'rec' - recenter back to input frame
		                  'ext' - cover both input and ref frame
	ext_pix             number of pixels to extend to save edge
	tmpdir              tmp file path
	------ OUTPUT ------
	'''
	def __init__(self, file, filREF=None, hdREF=None, \
		fmod='ref', ext_pix=0, tmpdir=None):
		'''
		self: hdr_ref, path_tmp, 
		(filIN, wmod, hdr, w, dim, Nx, Ny, Nw, im, wvl)
		'''
		## Set path of tmp files
		if tmpdir is None:
			path_tmp = os.getcwd()+'/tmp_proc/'
		else:
			path_tmp = tmpdir
		if not os.path.exists(path_tmp):
			os.makedirs(path_tmp)
		self.path_tmp = path_tmp

		## Inputs
		self.file = file
		self.filREF = filREF
		self.hdREF = hdREF
		self.fmod = fmod
		self.ext_pix = ext_pix
		
		## Init ref header
		self.hdr_ref = None

	def make_header(self, filIN, filREF=None, hdREF=None, fmod='ref', ext_pix=0):
		'''
		Make header tool

		------ INPUT ------
		filIN               single FITS file
		'''
		super().__init__(filIN)

		## Prepare reprojection header
		if filREF is not None:
			hdREF = read_fits(filREF).header
			# hdREF['EQUINOX'] = 2000.0

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
				w = fixwcs(header=hdREF).wcs
				try:
					pix_new = w.all_world2pix(world_arr, 1)
				except wcs.wcs.NoConvergence as e:
					pix_new = e.best_solution
					print("Best solution:\n{0}".format(e.best_solution))
					print("Achieved accuracy:\n{0}".format(e.accuracy))
					print("Number of iterations:\n{0}".format(e.niter))
				xmin = min(pix_new[:,0])
				xmax = max(pix_new[:,0])
				ymin = min(pix_new[:,1])
				ymax = max(pix_new[:,1])

				## Modify ref header
				if fmod=='rec': 
					hdREF['CRPIX1'] += -xmin
					hdREF['CRPIX2'] += -ymin
					hdREF['NAXIS1'] = math.ceil(xmax - xmin)
					hdREF['NAXIS2'] = math.ceil(ymax - ymin)
				elif fmod=='ext':
					if xmin<0:
						hdREF['CRPIX1'] += -xmin
					if ymin<0:
						hdREF['CRPIX2'] += -ymin
					hdREF['NAXIS1'] = math.ceil(max(xmax, hdREF['NAXIS1']-xmin, \
						xmax-xmin, hdREF['NAXIS1'])) + ext_pix # save edges
					hdREF['NAXIS2'] = math.ceil(max(ymax, hdREF['NAXIS2']-ymin, \
						ymax-ymin, hdREF['NAXIS2'])) + ext_pix
			## Save hdREF
			self.hdr_ref = hdREF

			## Test hdREF (Quick check: old=new or old<new)
			# w_new = fixwcs(header=hdREF).wcs
			# print('old: ', w.all_world2pix(
			# 	self.hdr['CRVAL1'], self.hdr['CRVAL2'], 1))
			# print('new: ', w_new.all_world2pix(
			# 	self.hdr['CRVAL1'], self.hdr['CRVAL2'], 1))
			# exit()
		else:
			print('ERROR: Can not find projection reference! ')
			exit()

	def make(self):
		'''
		Preparation (make header)
		'''
		file = self.file
		filREF = self.filREF
		hdREF = self.hdREF
		fmod = self.fmod
		ext_pix = self.ext_pix

		if isinstance(file, str):
			self.make_header(file, filREF, hdREF, fmod, ext_pix)
		elif isinstance(file, list):
			self.make_header(file[0], filREF, hdREF, fmod, ext_pix)
			if fmod=='ext':
				## Refresh self.hdr_ref in every circle
				for f in file:
					self.make_header(filIN=f, filREF=None, \
						hdREF=self.hdr_ref, fmod='ext', ext_pix=ext_pix)
		
		tqdm.write('<imontage> Making ref header...[done]')

		return self.hdr_ref

	def footprint(self, filOUT=None):
		'''
		Save reprojection footprint
		'''
		if filOUT is None:
			filOUT = self.path_tmp+'footprint'
		
		Nx = self.hdr_ref['NAXIS1']
		Ny = self.hdr_ref['NAXIS2']
		im_fp = np.ones((Ny, Nx))
		
		comment = "<imontage> footprint"
		write_fits(filOUT, self.hdr_ref, im_fp, COMMENT=comment)

		return im_fp

	def reproject(self, filIN, filOUT=None, \
		filUNC=None, dist='norm', postfix=''):
		'''
		Reproject 2D image or 3D cube

		------ INPUT ------
		filIN               single FITS file to reproject
		filOUT              output FITS file
		filUNC              unc files
		dist                uncertainty distribution
		                      'norm' - N(0,1)
		                      'splitnorm' - SN(0,lam,lam*tau)
		postfix              
		------ OUTPUT ------

		'''
		super().__init__(filIN)
		
		if dist=='norm':
			self.rand_norm(filUNC)
		elif dist=='splitnorm':
			self.rand_splitnorm(filUNC)
		
		## Set reprojection tmp path
		##---------------------------
		filename = os.path.basename(filIN)
		rep_tmp = self.path_tmp+filename+postfix+'/'
		if not os.path.exists(rep_tmp):
			os.makedirs(rep_tmp)

		self.slist = self.slice(rep_tmp+'slice', '_') # gauss_noise inclu

		## Do reprojection
		##-----------------
		cube_rep = []
		# for k in range(self.Nw):
			# hdr = self.hdr.copy()
			# for kw in self.hdr.keys():
			# 	if '3' in kw:
			# 		del hdr[kw]
			# hdr['NAXIS'] = 2
			# phdu = fits.PrimaryHDU(header=hdr, data=self.im[k,:,:])
			# im_rep = reproject_interp(phdu, self.hdr_ref)[0]
		for s in self.slist:
			im_rep = reproject_interp(s+'.fits', self.hdr_ref)[0]
			cube_rep.append(im_rep)
			write_fits(s+'rep_', self.hdr_ref, im_rep)
			fclean(s+'.fits')
		self.im = np.array(cube_rep)

		comment = "Reprojected by <imontage>. "
		if filOUT is None:
			filOUT = self.path_tmp+filename+postfix+'_rep'
		self.file_rep = filOUT
		
		write_fits(filOUT, self.hdr_ref, self.im, self.wvl, wmod=0, \
				COMMENT=comment)
		
		return self.im

	def reproject_mc(self, filIN, filUNC, Nmc=0, dist='norm'):
		'''
		Generate Monte-Carlo uncertainties for reprojected input file
		'''
		dataset = type('', (), {})()

		hyperim = [] # [j,(w,)y,x]
		for j in trange(Nmc+1, leave=False, \
			desc='<imontage> Reprojection (MC level)'):

			if j==0:
				im0 = self.reproject(filIN, \
					filUNC=None, dist=dist)
				file_rep = self.file_rep
			else:
				hyperim.append(self.reproject(filIN, \
					filUNC=filUNC, dist=dist, postfix='_'+str(j)))
		im0 = np.array(im0)
		hyperim = np.array(hyperim)
		unc = np.nanstd(hyperim, axis=0)

		comment = "An <imontage> production"
		
		write_fits(file_rep+'_unc', self.hdr_ref, unc, self.wvl, \
			COMMENT=comment)

		dataset.im0 = im0
		dataset.unc = unc
		dataset.hyperim = hyperim

		return dataset

	def combine(self, file, filOUT=None, method='avg', \
		filUNC=None, do_rep=True, Nmc=0, dist='norm'):
		'''
		Stitching input files (with the same wavelengths) to the ref WCS

		If filUNC is None, no MC
		If Nmc==0, no MC
		'''
		dataset = type('', (), {})()
		wvl = read_fits(file[0]).wave
		dataset.wvl = wvl

		superim0 = [] # [i,(w,)y,x]
		superunc = [] # [i,(w,)y,x]
		superim = [] # [i,j,(w,)y,x]
		Nf = np.size(file)
		for i in trange(Nf, leave=False, \
			desc='<imontage> Reprojection (file level)'):
			## (Re)do reprojection
			##---------------------
			if do_rep==True:
				## With MC
				if filUNC is not None:
					im0, unc, hyperim = self.reproject_mc(file[i], filUNC[i], \
						Nmc=Nmc, dist=dist)
					superunc.append(unc)
					superim.append(hyperim)
				## Without MC
				else:
					im0 = self.reproject(file[i])					
				superim0.append(im0)

			## Read archives
			##---------------
			else:
				filename = os.path.basename(file[i])
				file_rep = self.path_tmp+filename+'_rep'
				if filUNC is not None:
					hyperim = [] # [j,(w,)y,x]
					for j in range(Nmc+1):
						if j==0:
							superunc.append(read_fits(file_rep+'_unc').data)
						else:
							file_rep = self.path_tmp+filename+'_'+str(j)+'_rep'
							hyperim.append(read_fits(file_rep).data)
					hyperim = np.array(hyperim)
					superim.append(hyperim)
				superim0.append(read_fits(file_rep).data)

		superim0 = np.array(superim0)
		superunc = np.array(superunc)
		superim = np.array(superim)

		## Combine images
		##----------------
		## Think about using 'try - except'
		if filUNC is not None:
			inv_var = 1./superunc**2
			hyperim_comb = []
			for j in trange(Nmc+1, leave=False, \
				desc='<imontage> Stitching'):
				if j==0:
					if method=='avg':
						im0_comb = nanavg(superim0, axis=0)
					elif method=='wgt_avg':
						im0_comb = nanavg(superim0, axis=0, weights=inv_var)
				else:
					if method=='avg':
						hyperim_comb.append(nanavg(superim[:,j-1], axis=0))
					elif method=='wgt_avg':
						hyperim_comb.append(
							nanavg(superim[:,j-1], axis=0, weights=inv_var))
			hyperim_comb = np.array(hyperim_comb)
			unc_comb = np.nanstd(hyperim_comb)
		else:
			## If no unc, inverse variance weighted mean not available
			im0_comb = nanavg(superim0, axis=0)

		if filOUT is not None:
			comment = "An <imontage> production"

			write_fits(filOUT, self.hdr_ref, im0_comb, wvl, \
				COMMENT=comment)
			write_fits(filOUT+'_unc', self.hdr_ref, im0_comb, wvl, \
				COMMENT=comment)
		
		dataset.im0_comb = im0_comb
		dataset.unc_comb = unc_comb
		dataset.hyperim_comb = hyperim_comb
		dataset.superim0 = superim0
		dataset.superunc = superunc
		dataset.superim = superim

		tqdm.write('<imontage> Combining images...[done]')
		
		return dataset

	def clean(self, file=None):
		if file is not None:
			fclean(file)
		else:
			fclean(self.path_tmp)

class iswarp(impro):
	'''
	SWarp drop-in image montage toolkit
	i means <impro>-based
	In competetion with its fully Python-based twin <imontage>

	------ INPUT ------
	file                FITS files used to make header
	center              center of output image frame
	                      None - contains all input fields
	                      str('hh:mm:ss,dd:mm:ss') - manual input RA,DEC
	pixscale            pixel scale
	                      None - median of pixscale at center input frames
	                      float() - in arcseconds
	verbose             default: True
	tmpdir              tmp file path
	------ OUTPUT ------
	'''
	def __init__(self, file, center=None, pixscale=None, \
		verbose=True, tmpdir=None):
		'''
		self: path_tmp, verbose
		(filIN, wmod, hdr, w, dim, Nx, Ny, Nw, im, wvl)
		'''
		self.verbose = verbose
		## Set path of tmp files
		if tmpdir is None:
			path_tmp = os.getcwd()+'/tmp_swp/'
		else:
			path_tmp = tmpdir
		if not os.path.exists(path_tmp):
			os.makedirs(path_tmp)
		
		self.path_tmp = path_tmp

		## Images
		image_files = ' '
		for i in range(len(file)):
			image = read_fits(file[i]).data
			hdr = fixwcs(file[i]).header
			if len(image.shape)==3:
				## Extract 1st frame of the cube
				file[i] = path_tmp+os.path.basename(file[i])+'_ref'
				write_fits(file[i], hdr, image[0])
			
			image_files += file[i]+'.fits '

		## Create config file
		SP.call('swarp -d > swarp.cfg', shell=True, cwd=path_tmp)
		
		## Config param list
		swarp_opt = ' -c swarp.cfg -HEADER_ONLY Y -IMAGEOUT_NAME coadd.head '
		if center is not None:
			swarp_opt += ' -CENTER_TYPE MANUAL -CENTER '+center
		if pixscale is not None:
			swarp_opt += ' -PIXELSCALE_TYPE MANUAL -PIXEL_SCALE '+str(pixscale)
		if verbose==False:
			swarp_opt += ' -VERBOSE_TYPE QUIET '
		
		## Run SWarp
		SP.call('swarp '+swarp_opt+image_files, shell=True, cwd=path_tmp)

		## Save ref header 
		# self.refheader = read_fits(path_tmp+'coadd').header
		self.refheader = None

	def footprint(self, filOUT=None):
		'''
		Save reprojection footprint
		'''
		if filOUT is None:
			filOUT = self.path_tmp+'footprint'
		
		Nx = self.refheader['NAXIS1']
		Ny = self.refheader['NAXIS2']
		im_fp = np.ones((Ny, Nx))
		
		comment = "<iswarp> footprint"
		write_fits(filOUT, self.refheader, im_fp, COMMENT=comment)

		return im_fp

	def combine(self, file, combtype='med', \
		keepedge=False, uncpdf=None, filOUT=None):
		'''
		Combine 

		------ INPUT ------
		filIN              input FITS files should have the same wvl
		combtype            combine type
		                      med - median
		                      avg - average
		                      wgt_avg - inverse variance weighted average
		keepedge            default: False
		filOUT              output FITS file
		------ OUTPUT ------
		'''
		cl = type('', (), {})()
		verbose = self.verbose
		path_tmp = self.path_tmp

		## Images and weights
		##--------------------
		Nf = len(file)
		
		imshape = read_fits(file[0]).data.shape
		if len(imshape)==3:
			Nw = imshape[0]
			wvl = read_fits(file[0]).wave
		else:
			Nw = 1
		
		## Build imlist & wgtlist (size=Nf)
		imlist = []
		wgtlist = []
		for i in range(Nf):
			filename = os.path.basename(file[i])
			## Set slice file
			path_slice = path_tmp+'slices/'
			if not os.path.exists(path_slice):
				os.makedirs(path_slice)
			file_slice = path_slice+filename
			## Slice
			super().__init__(file[i])
			if uncpdf=='norm':
				self.rand_norm(file[i]+'_unc')
			elif uncpdf=='splitnorm':
				self.rand_splitnorm(file[i]+'_unc')
			imlist.append(self.slice(file_slice, ''))
			
			if combtype=='wgt_avg':
				super().__init__(file[i]+'_unc')
				wgtlist.append(self.slice_inv_sqrt(file_slice, '.weight'))

		## Build image_files & weight_files (size=Nw)
		image_files = [' ']*Nw
		weight_files = [' ']*Nw

		## Let's SWarp
		##-------------
		hyperimage = []
		for k in trange(Nw, leave=False, 
			desc='<iswarp> Combining (by wvl)'):
			for i in range(Nf):
				image_files[k] += imlist[i][k]+'.fits '

				if combtype=='wgt_avg':
					weight_files[k] += wgtlist[i][k]+'.fits '

			## Create config file
			SP.call('swarp -d > swarp.cfg', shell=True, cwd=path_tmp)
			## Config param list
			swarp_opt = ' -c swarp.cfg -SUBTRACT_BACK N '
			if combtype=='med':
				pass
			elif combtype=='avg':
				swarp_opt += ' -COMBINE_TYPE AVERAGE '
			elif combtype=='wgt_avg':
				swarp_opt += ' -COMBINE_TYPE WEIGHTED '
				swarp_opt += ' -WEIGHT_TYPE MAP_WEIGHT '
				swarp_opt += ' -WEIGHT_SUFFIX .weight.fits '
				# swarp_opt += ' -WEIGHT_IMAGE '+weight_files[k] # not worked
			if verbose==False:
				swarp_opt += ' -VERBOSE_TYPE QUIET '
			## Run SWarp
			SP.call('swarp '+swarp_opt+' -RESAMPLING_TYPE LANCZOS3 '+image_files[k], \
				shell=True, cwd=path_tmp)
			coadd = read_fits(path_tmp+'coadd')
			newimage = coadd.data
			newheader = coadd.header

			## Add back in the edges because LANCZOS3 kills the edges
			## Do it in steps of less and less precision
			if keepedge==True:
				oldweight = read_fits(path_tmp+'coadd.weight').data
				if np.sum(oldweight==0)!=0:
					SP.call('swarp '+swarp_opt+' -RESAMPLING_TYPE LANCZOS2 '+image_files[k], \
						shell=True, cwd=path_tmp)
					edgeimage = read_fits(path_tmp+'coadd').data
					newweight = read_fits(path_tmp+'coadd.weight').data
					edgeidx = np.ma.array(oldweight, 
						mask=np.logical_and(oldweight==0, newweight!=0)).mask
					if edgeidx.any():
						newimage[edgeidx] = edgeimage[edgeidx]

					oldweight = read_fits(path_tmp+'coadd.weight').data
					if np.sum(oldweight==0)!=0:
						SP.call('swarp '+swarp_opt+' -RESAMPLING_TYPE BILINEAR '+image_files[k], \
							shell=True, cwd=path_tmp)
						edgeimage = read_fits(path_tmp+'coadd').data
						newweight = read_fits(path_tmp+'coadd.weight').data
						edgeidx = np.ma.array(oldweight, 
							mask=np.logical_and(oldweight==0, newweight!=0)).mask
						if edgeidx.any():
							newimage[edgeidx] = edgeimage[edgeidx]

						oldweight = read_fits(path_tmp+'coadd.weight').data
						if np.sum(oldweight==0)!=0:
							SP.call('swarp '+swarp_opt+' -RESAMPLING_TYPE NEAREST '+image_files[k], \
								shell=True, cwd=path_tmp)
							edgeimage = read_fits(path_tmp+'coadd').data
							newweight = read_fits(path_tmp+'coadd.weight').data
							edgeidx = np.ma.array(oldweight, 
								mask=np.logical_and(oldweight==0, newweight!=0)).mask
							if edgeidx.any():
								newimage[edgeidx] = edgeimage[edgeidx]
			
			## Astrometric flux-rescaling based on the local ratio of pixel scale
			## Complementary for lack of FITS kw 'FLXSCALE'
			## Because SWarp is conserving surface brightness/pixel
			oldcdelt = pc2cd(wcs=fixwcs(file[i]).wcs).cdelt
			newcdelt = [newheader['CD1_1'], newheader['CD2_2']]
			old_pixel_fov = abs(oldcdelt[0]*oldcdelt[1])
			new_pixel_fov = abs(newcdelt[0]*newcdelt[1])
			newimage = newimage * old_pixel_fov/new_pixel_fov
			zeroma = np.ma.array(newimage, mask=(newimage==0)).mask
			newimage[zeroma] = np.nan
			write_fits(path_slice+'coadd_'+str(k), newheader, newimage)
			hyperimage.append(newimage)

		hyperimage = np.array(hyperimage)
		if filOUT is not None:
			write_fits(filOUT, newheader, hyperimage, wvl)

		cl.image = hyperimage 
		cl.header = newheader

		return cl

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
	filOUT              output file
	------ OUTPUT ------
	'''
	def __init__(self, filIN, filKER, kfile, \
		filUNC=None, dist='norm', psf=None, filTMP=None, filOUT=None):
		## INPUTS
		super().__init__(filIN)
		
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
			write_fits(self.filOUT, self.hdr, self.convim, self.wvl, \
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
	parobs[2]           IRC N3 (long exp) frame (2MASS corrected; 90 deg rot)
	Nw                  num of wave
	Ny                  slit length
	Nx                  slit width
	------ OUTPUT ------
	'''
	def __init__(self, ipath=None, parobs=None, verbose=False):
		self.path = ipath + parobs[0] + '/irc_specred_out_' + parobs[1]+'/'
		filIN = self.path + parobs[2]
		super().__init__(filIN)
		
		self.filSAV = self.path + parobs[0] + '.N3_NG.IRC_SPECRED_OUT'
		self.table = readsav(self.filSAV+savext, python_dict=True)['source_table']

		## Slit width will be corrected by intercalib with IRS data
		if parobs[1]=='Ns':
			self.slit_width = 3 # 412pix * 5"/10' = 3.43 pix (Ns)
		elif parobs[1]=='Nh':
			self.slit_width = 2 # 412pix * 3"/10' = 2.06 pix (Nh)

		if verbose==True:
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

	def spec_build(self, filOUT=None, write_unc=True, Ny=31, sig_pt=0.):
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

			if write_unc==True:
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
	oldimage = ds.data
	hdr = ds.header
	w = fixwcs(filIN).wcs
	# hdr['NAXIS1'] = x1 - x0 + 1
	# hdr['NAXIS2'] = y1 - y0 + 1
	hdr['CRPIX1'] += -x0
	hdr['CRPIX2'] += -y0
	newimage = oldimage[y0:y1+1, x0:x1+1]

	write_fits(filOUT, hdr, newim)

	return newimage

def hswarp(oldimage, oldheader, refheader, \
	keepedge=False, tmpdir=None, verbose=True):
	'''
	Python version of hswarp (IDL), 
	a SWarp drop-in replacement for hastrom, 
	created by S. Hony

	------ INPUT ------
	oldimage            ndarray
	oldheader           header object
	refheader           ref header
	keepedge            default: False
	tmpdir              default: None
	verbose             default: True
	------ OUTPUT ------
	cl                  output object
	  image               newimage
	  header              newheader
	'''
	## Initialize output object
	cl = type('', (), {})()

	## Set path of tmp files
	if tmpdir is None:
		path_tmp = os.getcwd()+'/tmp_hswarp/'
	else:
		path_tmp = tmpdir
	if not os.path.exists(path_tmp):
		os.makedirs(path_tmp)

	## Make input
	write_fits(path_tmp+'old', oldheader, oldimage)
	with open(path_tmp+'coadd.head', 'w') as f:
		f.write(str(refheader))

	## Create config file
	SP.call('swarp -d > swarp.cfg', shell=True, cwd=path_tmp)
	## Config param list
	swarp_opt = ' -c swarp.cfg -SUBTRACT_BACK N '
	if verbose==False:
		swarp_opt += ' -VERBOSE_TYPE QUIET '
	## Run SWarp
	SP.call('swarp '+swarp_opt+' -RESAMPLING_TYPE LANCZOS3 '+' old.fits', \
		shell=True, cwd=path_tmp)
	coadd = read_fits(path_tmp+'coadd')
	newimage = coadd.data
	newheader = coadd.header

	## Add back in the edges because LANCZOS3 kills the edges
	## Do it in steps of less and less precision
	if keepedge==True:
		oldweight = read_fits(path_tmp+'coadd.weight').data
		if np.sum(oldweight==0)!=0:
			SP.call('swarp '+swarp_opt+' -RESAMPLING_TYPE LANCZOS2 '+' old.fits', \
				shell=True, cwd=path_tmp)
			edgeimage = read_fits(path_tmp+'coadd').data
			newweight = read_fits(path_tmp+'coadd.weight').data
			edgeidx = np.ma.array(oldweight, 
				mask=np.logical_and(oldweight==0, newweight!=0)).mask
			if edgeidx.any():
				newimage[edgeidx] = edgeimage[edgeidx]

			oldweight = read_fits(path_tmp+'coadd.weight').data
			if np.sum(oldweight==0)!=0:
				SP.call('swarp '+swarp_opt+' -RESAMPLING_TYPE BILINEAR '+' old.fits', \
					shell=True, cwd=path_tmp)
				edgeimage = read_fits(path_tmp+'coadd').data
				newweight = read_fits(path_tmp+'coadd.weight').data
				edgeidx = np.ma.array(oldweight, 
					mask=np.logical_and(oldweight==0, newweight!=0)).mask
				if edgeidx.any():
					newimage[edgeidx] = edgeimage[edgeidx]

				oldweight = read_fits(path_tmp+'coadd.weight').data
				if np.sum(oldweight==0)!=0:
					SP.call('swarp '+swarp_opt+' -RESAMPLING_TYPE NEAREST '+' old.fits', \
						shell=True, cwd=path_tmp)
					edgeimage = read_fits(path_tmp+'coadd').data
					newweight = read_fits(path_tmp+'coadd.weight').data
					edgeidx = np.ma.array(oldweight, 
						mask=np.logical_and(oldweight==0, newweight!=0)).mask
					if edgeidx.any():
						newimage[edgeidx] = edgeimage[edgeidx]

	## SWarp is conserving surface brightness/pixel
	## while the pixels size changes
	oldcdelt = pc2cd(wcs=fixwcs(header=oldheader).wcs).cdelt
	refcdelt = pc2cd(wcs=fixwcs(header=refheader).wcs).cdelt
	old_pixel_fov = abs(oldcdelt[0]*oldcdelt[1])
	new_pixel_fov = abs(refcdelt[0]*refcdelt[1])
	newimage = newimage * old_pixel_fov/new_pixel_fov
	print('-------------------')
	print(old_pixel_fov/new_pixel_fov)
	write_fits(path_tmp+'new', newheader, newimage)
	print('-------------------')

	## Delete tmp file if tmpdir not given
	if tmpdir is None:
		fclean(path_tmp)

	cl.image = newimage
	cl.header = newheader

	return cl

def concatenate(flist, filOUT=None, comment=None, sort_wave=True):
	'''

	'''
	dataset = type('', (), {})()

	## Read data
	wave = []
	data = []
	for f in flist:
		ds = read_fits(f)
		data.append(ds.data)
		wave.append(ds.wave)
	data = np.concatenate(data)
	wave = np.concatenate(wave)
	hdr = ds.header

	## Sort wavelength in wave ascending order
	if sort_wave==True:
		ind = sorted(range(len(wave)), key=wave.__getitem__)
		wave = np.sort(wave)
		data = data[ind]

	## Save data
	dataset.data = data
	dataset.wave = wave

	## Write FITS file
	if filOUT is not None:
		write_fits(filOUT, hdr, data, wave, \
			COMMENT=comment)

	return dataset

"""
------------------------------ MAIN (test) ------------------------------
"""
if __name__ == "__main__":

	import os
	path_test = os.getcwd() + '/../tests/'
