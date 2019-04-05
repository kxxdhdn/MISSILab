#!/usr/bin/env python
# -*- coding: utf-8 -*-

from astropy.io import fits
from astropy import units as u
from astropy import wcs
from astropy.modeling import models, fitting
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import Axes3D

"""
------------------------------ FUNCTION ------------------------------
"""
def bimapview(lx, ly, wvl, data, data0, wvlnum, patch):
	"""
	view the map and choose 1 point/region to check spectrum & flux distribution

	--- INPUT ---
	lx           length of map
	ly           width of map
	wvl          wavelength
	data         cube
	data0        reference cube
	wvlnum       number of cube slice, int type
	patch        replot in order to add selected patch
	--- OUTPUT ---
	fig1         current figure
	"""

	x = np.arange(lx)
	y = np.arange(ly)
	fig1 = plt.figure("map", figsize=(14,5))
	fig1.suptitle("IRS maps at {} micron".format(wvl[wvlnum]))
	cnorm = None # Linear
#	cnorm = colors.PowerNorm(gamma=0.5) # Square root
#	cnorm = colors.SymLogNorm(linthresh=1., linscale=1., clip=False) # Logarithmic
	# ----- target cube -----
	ax11 = plt.subplot(121)
	pcm = ax11.pcolormesh(x, y, data[wvlnum,:,:], norm=cnorm, cmap='rainbow')
	fig1.colorbar(pcm, ax=ax11, extend='max')
	ax11.set_title("dh")
	ax11.set_xlabel("x")
	ax11.set_ylabel("y")
	# ----- ref cube -----
	ax12 = plt.subplot(122)
	pcm = ax12.pcolormesh(x, y, data0[wvlnum,:,:], norm=cnorm, cmap='rainbow')
	fig1.colorbar(pcm, ax=ax12, extend='max')
	ax12.set_title("sh")
	ax12.set_xlabel("x")
	ax12.set_ylabel("y")
	if np.size(patch)==4:
		## add rectangle patch
		plt.gca().add_patch(plt.Rectangle(
			(patch[0], patch[1]), # (xmin, ymin)
			patch[2]-patch[0], patch[3]-patch[1], # width, height
    		edgecolor='red', linestyle='-', fill=False))
	elif np.size(patch)==3:
		## add circle patch
		plt.gca().add_patch(plt.Circle(
			(patch[0], patch[1]), # (x0, y0)
			patch[2], # radius
    		edgecolor='red', linestyle='-', fill=False))

	return fig1

def bispecview(wvl, data1, data2):
	"""
	extract spectrum at chosen square

	--- INPUT ---
	wvl           wavelength
	data1         box
	data2         box in ref frame
	"""
	fig2 = plt.figure("box spectrum")
	ax2 = plt.subplot()
	ax2.errorbar(wvl, np.nansum(data1, axis=0), 0.,
		c='k', ls='-', lw=.5, ecolor='r', elinewidth=1., label="data")
	ax2.errorbar(wvl, np.nansum(data2, axis=0), 0.,
		c='b', ls='-', lw=.5, ecolor='r', elinewidth=1., label="ref")
	ax2.set_title("Spectra")
	ax2.set_xlabel('Wavelength (micron)')
	ax2.set_ylabel(r'$F_{\nu} (MJy/sr)$')
#	ax2.set_yscale('log', nonposy='clip')
	ax2.legend(loc='upper left')
#	fig2.savefig('spectrum.png')

def cal_res(wvl, data1, data2, trim):
	"""
	calculate residuals for chosen pixels with each wavelength

	--- INPUT ---
	wvl       
	data1        data (2D)
	data2        ref
	trim         alert trim
	--- OUTPUT ---
	res          1D array
	res_std      std_dev of res
	npix         nb of pixels used in the calculation
	badv         wvl of slice which might have bad pixels ()
	"""
	npix = np.size(data1[:,0])
	print("npix = ", npix)
	nv = np.size(wvl)
	res = np.zeros(nv, dtype=float)
	badv = []
	for k in range(nv):
		F_dh = np.nansum(data1[:,k])
		F_sh = np.nansum(data2[:,k])
		if F_sh==0.:
			res[k] = 0.
		else:
			res[k] = (F_dh - F_sh) / F_sh
	## alert
		if abs(res[k])>trim:
			badv.append(wvl[k])
			print("Need to review {:.5} micron slice".format(wvl[k]))
	print(badv)
	## res_std (par rapport a wvl)
	res_std = rms(res)
	## plot res - wvl
	fig3 = plt.figure("Residual Distribution")
	ax3 = plt.subplot()
	ax3.plot(wvl, res, 'k-', lw=.5, label="std_dev = {:.2%}".format(res_std))
	ax3.set_xlabel('Wavelength (micron)')
	ax3.set_ylabel('Residual')
#	ax3.set_xscale('log', nonposy='clip')
	ax3.legend(loc='upper right')

	return res

def rms(a):
	"""
	calculate root mean square
	
	--- INPUT ---
	a            an array
	--- OUTPUT ---
	rms          root mean square of a
	"""
	ms = np.nanmean(a*a)

	return np.sqrt(ms)


"""
------------------------------ MAIN ------------------------------
"""
filename = '1_LL2_cube'
filename0 = 'n83_84_LL2_cube' 
#filename0 = filename # set 2 filenames the same if no need to compare cubes
## select fixed rectangle: lower left (x1, y1) & upper right (x2, y2)
#xmin, xmax = 50, 70
#ymin, ymax = 65, 80
xmin, xmax = 52, 54
ymin, ymax = 67, 69

with fits.open(filename+'.fits') as hdul:
	data = hdul[0].data
	hdr = hdul[0].header
	NAXIS1 = hdr['NAXIS1']
	NAXIS2 = hdr['NAXIS2']
	NAXIS3 = hdr['NAXIS3']
	wvl = hdul[1].data[0][0][:,0]
#	wvl = hdul[1].data # for rewitten header
hdul0 = fits.open(filename0+'.fits')
data0 = hdul0[0].data

while True:
	b1 = input("(Re)Analyse? (0 - yes / else - no) ") # keyboard input 1
	if b1=="0":
		wvlnum = int(input("Select a slice of wavelength (0 - {}): ".format(NAXIS3-1))) # keyboard input 2
		print("wavelength = ", wvl[wvlnum])

		## view map 
		figmap = bimapview(NAXIS1, NAXIS2, wvl, data, data0, wvlnum, 0)
		"""
		box statistics
		"""
		## select box mode
		b0 = input("Select mode (0 - fixed / 1 - clicked / 2 - circle / else - skip) ") # keyboard input 0
		print("...")
		if b0=="0":
			## fixed box 
			xb, xt = xmin, xmax # b-bottom, t-top
			yb, yt = ymin, ymax
			plt.close(figmap)
		elif b0=="1":
			## left click to select 2 points which determine a rectangle by the diagonal
			pts = figmap.ginput(n=2, show_clicks=True)
			(x1, y1) = np.array(pts[0]).astype(int)
			(x2, y2) = np.array(pts[1]).astype(int)
			xb, xt = min(x1, x2), max(x1, x2)
			yb, yt = min(y1, y2), max(y1, y2)
			plt.close(figmap)
		elif b0=="2":
			## left click to select 1 point which determine a circle
			radius = int(input("Select circle radius (only int type): ")) # keyboard input 3
			pts = figmap.ginput(n=1, show_clicks=True)
			(x0, y0) = np.array(pts[0]).astype(int)
			plt.close(figmap)
		else:
			plt.close(figmap)
			continue
		
		plt.show()

		## extract the flux
		if not(b0=="2"):
			print("rectangle diagonal ({}, {}) - ({}, {})".format(xb, yb, xt, yt))
			## add patch onto map
			patch = [xb, yb, xt, yt]
			bimapview(NAXIS1, NAXIS2, wvl, data, data0, wvlnum, patch)
			## fluxes in box 
			data1 = data[:, yb:yt, xb:xt].reshape(((xt-xb)*(yt-yb), NAXIS3))
			data2 = data0[:, yb:yt, xb:xt].reshape(((xt-xb)*(yt-yb), NAXIS3))
		elif b0=="2":
			print("circle centered at ({0}, {1})".format(x0, y0))
			## add patch onto map
			patch = [x0, y0, radius]
			bimapview(NAXIS1, NAXIS2, wvl, data, data0, wvlnum, patch)
			## fluxes in box 
			data1 = [data[:, y0, x0]]
			data2 = [data0[:, y0, x0]]
			for dx in range(radius+1):
				for dy in range(radius+1):
					p2 = dx*dx + dy*dy
					r2 = radius*radius
					if (dx==0 and dy==0):
						continue
					elif (dx==0 or dy==0):
						#print("dx, dy ", dx, dy)
						data1.append(data[:, y0-dy, x0-dx])
						data1.append(data[:, y0+dy, x0+dx])
						data2.append(data0[:, y0-dy, x0-dx])
						data2.append(data0[:, y0+dy, x0+dx])
					elif (p2!=0 and p2<=r2):
						#print("dx, dy ", dx, dy)
						data1.append(data[:, y0-dy, x0-dx])
						data1.append(data[:, y0-dy, x0+dx])
						data1.append(data[:, y0+dy, x0-dx])
						data1.append(data[:, y0+dy, x0+dx])
						data2.append(data0[:, y0-dy, x0-dx])
						data2.append(data0[:, y0-dy, x0+dx])
						data2.append(data0[:, y0+dy, x0-dx])
						data2.append(data0[:, y0+dy, x0+dx])
			data1 = np.array(data1)
			data2 = np.array(data2)
			print("(nb of pix in circle, NAXIS3) = ", data1.shape)

		## view spectrum
		bispecview(wvl, data1, data2)
		## compare two cubes [comment the following if not compare]
		res = cal_res(wvl, data1, data2, .5)

		plt.show()
	else:
		exit(0)
