class imontage(impro):
	'''
	2D image or 3D cube montage toolkit
	i means <impro>-based or initialize

	------ INPUT ------
	file                input FITS file (list)
	filREF              ref file (priority if co-exist with input header)
	hdREF               ref header
	fmod                output file frame
		                  'ref' - same as ref frame (Default)
		                  'rec' - recenter back to input frame
		                  'ext' - cover both input and ref frame
	ext_pix             number of pixels to extend to save edge
	------ OUTPUT ------
	'''
	def __init__(self, file, filREF=None, hdREF=None, fmod='ref', ext_pix=0):
		'''
		self: hdr_ref, path_tmp, 
		(filIN, wmod, hdr, w, dim, Nx, Ny, Nw, im, wvl)
		'''
		## Set path of tmp files
		path_tmp = os.getcwd()+'/tmp_proc/'
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

	def make_header(self, file, filREF=None, hdREF=None, fmod='ref', ext_pix=0):
		'''
		Make header tool

		------ INPUT ------
		file                single FITS file
		'''
		super().__init__(file)

		## Prepare reprojection header
		if filREF is not None:
			hdREF = ext_wcs(ref).header
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
						xmax-xmin, hdREF['NAXIS1'])) + ext_pix # save edges
					hdREF['NAXIS2'] = math.ceil(max(ymax, hdREF['NAXIS2']-ymin, \
						ymax-ymin, hdREF['NAXIS2'])) + ext_pix
			## Save hdREF
			self.hdr_ref = hdREF
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
					self.make_header(file=f, filREF=None, \
						hdREF=self.hdr_ref, fmod='ext', ext_pix=ext_pix)

		return self.hdr_ref

	def footprint(self, filOUT=None, wmod=0):
		'''
		Show reprojection footprint
		'''
		if filOUT is None:
			filOUT = self.path_tmp+'footprint'
		
		Nx = self.hdr_ref['NAXIS1']
		Ny = self.hdr_ref['NAXIS2']
		im_fp = np.ones((Ny, Nx))
		
		comment = "<imontage> footprint"
		write_fits(filOUT, self.hdr_ref, im_fp, None, wmod, \
				COMMENT=comment)

		return im_fp

	def reproject(self, filOUT=None, filUNC=None, dist='norm', wmod=0, filTMP=None):
		'''
		Reproject 2D image or 3D cube

		------ INPUT ------
		filOUT              output FITS file
		wmod                wave mode of output file
		filUNC              unc files
		dist                uncertainty distribution
		                      'norm' - N(0,1)
		                      'splitnorm' - SN(0,lam,lam*tau)
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
		if filTMP is None:
			fclean(self.path_tmp)

		self.im = np.array(cube_rep)

		if filOUT is not None:
			self.hdr = self.hdr_ref

			comment = "Reprojected by <iproject>. "
			write_fits(filOUT, self.hdr, self.im, self.wvl, self.wmod, \
				COMMENT=comment)
	
		return self.im

	def reproject_mc(self, Nmc=0, filUNC=None, dist='norm', write_mc=False):
		'''
		Generate Monte-Carlo uncertainties
		'''
		file_rep = None
		hyperim = [] # [j,(w,)y,x]
		for j in range(Nmc+1):
			if write_mc==True:
				if not os.path.exists(f+'/'):
					os.makedirs(f+'/')
				file_rep = f + '/rep_' + str(j)
				file_unc = f + '/unc'

			if j==0:
				im = self.reproject(filOUT=file_rep, \
					filUNC=None, dist=dist, wmod=0, filTMP=None)
			else:
				hyperim.append(self.reproject(filOUT=file_rep, \
					filUNC=filUNC, dist=dist, wmod=0, filTMP=None))
		hyperim = np.array(hyperim)
		sigma = np.nanstd(hyperim, axis=0)
			
		return im, sigma

	def combine(self, filOUT=None, method='average', write_mc=False, \
		do_rep=True, Nmc=0, ulist=None, dist='norm'):
		'''
		Stitching input files (with the same wavelengths) to the ref WCS

		'''
		hyperim = [] # [i,(w,)y,x]
		unc = [] # [i,(w,)y,x]]
		for i, f in enumerate(self.file):
			if do_rep==True:
				im, sigma = self.reproject_mc(Nmc=Nmc, \
					filUNC=ulist[i], dist=dist, write_mc=write_mc)
				hyperim.append(im)
				unc.append(sigma)
			else:
				ds = read_fits(f)
				hyperim.append(ds.data)
				if ulist is not None:
					unc.append(read_fits(ulist[i]).data)
		hyperim = np.array(hyperim)
		unc = np.array(unc)

		## Combine images
		if method=='average':
			im_comb = nanavg(hyperim, axis=0)
		elif method=='weighted_avg':
			im_comb = nanavg(hyperim, axis=0, weights=unc)

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
