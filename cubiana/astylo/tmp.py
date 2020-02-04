class iswarp(impro):
	'''
	SWarp drop-in image montage toolkit
	i means <impro>-based
	In competetion with its fully Python-based twin <imontage>

	------ INPUT ------
	file                FITS file
	center              center of output image frame
	                      None - contains all input fields
	                      str('hh:mm:ss,dd:mm:ss') - manual input RA,DEC
	pixscale            pixel scale
	                      'med' - median of pixscale at center input frames
	                      float() - in arcseconds
	verbose             default: True
	tmpdir              tmp file path
	------ OUTPUT ------
	'''
	def __init__(self, file, center=None, pixscale='med', \
		verbose=True, tmpdir=None):
		'''
		self: file, refheader, path_tmp, verbose
		(filIN, wmod, hdr, w, dim, Nx, Ny, Nw, im, wvl)
		'''
		self.file = file
		self.verbose = verbose
		## Set path of tmp files
		if tmpdir is None:
			path_tmp = os.getcwd()+'/tmp_proc/'
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
		if pixscale=='med':
			pass
		else:
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

	def reproject(self, , combtype='med', keepedge=False, ):
		'''
		reproject

		------ INPUT ------
		filOUT              output FITS file
		combtype            combine type
		                      med - median
		                      avg - average
		                      wgt_avg - inverse variance weighted average
		keepedge            default: False
		------ OUTPUT ------
		'''
		cl = type('', (), {})()

		image_files += f+'.fits '
		weight_files += f+'_wgt.fits '

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
			swarp_opt += ' -WEIGHT_TYPE MAP_WEIGHT'
			swarp_opt += ' -WEIGHT_IMAGE '+weight_files
		if verbose==False:
			swarp_opt += ' -VERBOSE_TYPE QUIET '
		
		## Run SWarp
		SP.call('swarp '+swarp_opt+' -RESAMPLING_TYPE LANCZOS3 '+image_files, \
			shell=True, cwd=path_tmp)
		coadd = read_fits()
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
		
		## Astrometric flux-rescaling based on the local ratio of pixel scale
		## Complementary for lack of FITS kw 'FLXSCALE'
		## Because SWarp is conserving surface brightness/pixel
		oldcdelt = pc2cd(wcs=fixwcs(header=oldheader).wcs).cdelt
		refcdelt = pc2cd(wcs=fixwcs(header=refheader).wcs).cdelt
		old_pixel_fov = abs(oldcdelt[0]*oldcdelt[1])
		new_pixel_fov = abs(refcdelt[0]*refcdelt[1])
		newimage = newimage * old_pixel_fov/new_pixel_fov
		# write_fits(path_tmp+'new', newheader, newimage)

		cl.image = newimage
		cl.header = newheader

		return cl

	def combine(self, filOUT):
		'''
		Combine

		'''
		cl = type('', (), {})()
		file = self.file
		verbose = self.verbose
		path_tmp = self.path_tmp

		## Images and weights
		image_files = ' '
		weight_files = ' '
		for i in trange(len(file), leave=False, \
			desc='<iswarp> Reading input files'):

			image = read_fits(file[i]).data
			hdr = fixwcs(file[i]).header
			if len(image.shape)==3:
				## Extract 1st frame of the cube
				file[i] = path_tmp+os.path.basename(file[i])+'_ref'
				write_fits(file[i], hdr, image[0])
			
			image_files += file[i]+'.fits '

		


def hswarp_pro(flist, refheader, combtype='med', \
	keepedge=False, tmpdir=None, verbose=True):

	## Initialize output object
	cl = type('', (), {})()

	## Delete tmp file if tmpdir not given
	if tmpdir is None:
		fclean(path_tmp)

	cl.image = newimage
	cl.header = newheader

	return cl
