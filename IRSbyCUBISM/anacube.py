#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
anacube.py is a simple routine dedicated to extracting spectra from cubes built by CUBISM 
for further analysis such as cube quality estimation by inspecting spatial flux distribution. 
It has a twin version anacube_dual.py where the comparison of two cubes are possible.


Comment fluxdist function, line ?, if you only want to view spectrum.
Change fluxdist input "wvlnum" where "all_wvl" option means taking all wavelengths into acount.

Modify filename and parameters at the beginning of MAIN procedure.

Uncertainty is not yet available to be added to spectrum.
"""

import os, sys
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
def mapview(lx, ly, wvl, data, wvlnum, patch):
	"""
	view the map and select 1 point/region to check spectrum & flux distribution

	--- INPUT ---
	lx           length of map
	ly           width of map
	wvl          wavelength
	data         cube
	wvlnum       number of cube slice, int type
	patch        replot in order to add selected patch
	--- OUTPUT ---
	fig1         current figure
	"""
	fig1 = plt.figure("map", figsize=(10,8))
	cnorm = None # Linear
#	cnorm = colors.PowerNorm(gamma=0.5) # Square root
#	cnorm = colors.SymLogNorm(linthresh=1., linscale=1., clip=False) # Logarithmic
	ax1 = plt.subplot()
	x = np.arange(lx)
	y = np.arange(ly)
	pcm = ax1.pcolormesh(x, y, data[wvlnum,:,:], norm=cnorm, cmap='rainbow')
	fig1.colorbar(pcm, ax=ax1, extend='max')
	ax1.set_title("IRS map at {} micron".format(wvl[wvlnum]))
	ax1.set_xlabel("x")
	ax1.set_ylabel("y")
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

def specview(wvl, data, unc):
	"""
	extract spectrum at chosen square

	--- INPUT ---
	wvl          wavelength
	data         box
	unc          recalculated uncertainty box
	"""
	fig2 = plt.figure("box spectrum")
	ax2 = plt.subplot()
	ax2.errorbar(wvl, np.nansum(data, axis=(1,2)), unc,
		c='k', ls='-', lw=.5, ecolor='r', elinewidth=1., label="IRS data")
	ax2.set_title("Spectrum")
	ax2.set_xlabel('Wavelength (micron)')
	ax2.set_ylabel(r'$F_{\nu} (MJy)$')
#	ax2.set_yscale('log', nonposy='clip')
	ax2.legend(loc='upper left')
#	fig2.savefig('spectrum.png')

def fluxdist(wvl, data, wvlnum):
	"""
	plot flux distribution for pixels in the box

	--- INPUT ---
	wvl          wavelength
	data         square
	"""
	fig3 = plt.figure("box flux distribution")
	ax3 = plt.subplot()
	if wvlnum=="all_wvl":
		data1 = data[:,:,:].reshape(-1)
	else:
		data1 = data[wvlnum,:,:].reshape(-1)
	num, bins, patches = ax3.hist(data1, bins=int(np.size(data1)),
		alpha=.3, lw=.5, edgecolor='w', color='c')
	bins_c = []
	for i in range(np.size(bins)-1):
		bins_c.append((bins[i]+bins[i+1])/2.)
	bins_c = np.array(bins_c)

	## Fit the data using a Gaussian
	g_init = models.Gaussian1D(amplitude=1., mean=0, stddev=1.)
	fit_g = fitting.LevMarLSQFitter()
	g = fit_g(g_init, bins_c, num)
	numsum = 0
	for j in range(np.size(num)):
		if (bins_c[j]>median-g.stddev and bins_c[j]<g.mean+g.stddev):
			numsum += num[j]
	sig1 = float(numsum / np.nansum(num))
	print("box flux mean: ", g.mean*1.)
	print("box flux median: ", )
	print("5-sigma from the median percentage: {:.0%}".format(sig1))
	ax3.plot(bins_c, g(bins_c), c='r')
	ax3.axvline(g.mean, c='k', ls='-.', 
		label="mean: {:.2}".format(g.mean*1.))
	ax3.axvline(g.mean-g.stddev, c='b', ls=':', 
		label="1-sigma ratio: {:.0%}".format(sig1))
	ax3.axvline(g.mean+g.stddev, c='b', ls=':')

	if wvlnum=="all_wvl":
		ax3.set_title("Flux distribution")
	else:
		ax3.set_title("Flux distribution at {} micron".format(wvl[wvlnum]))
	ax3.set_xlabel(r'$F_{\nu} (MJy)$')
	ax3.set_ylabel("Number of pixels")
	ax3.legend(loc='upper right')

def calunc(wvl, data, unc):
	"""
	recalculate uncertainties for pixels in the box

	--- INPUT ---
	wvl          wavelength
	data         box
	unc          uncertainty box (3D)
	--- OUTPUT ---
	unc_wvl      1D unc
	"""
	unc_wvl = unc

	return unc_wvl

"""
------------------------------ MAIN ------------------------------
"""
filename = 'n66north_SL2_cube'
uncname = 'n66north_SL2_cube_unc'
## select fixed rectangle: lower left (x1, y1) & upper right (x2, y2)
xmin, xmax = 90, 120
ymin, ymax = 70, 100

with fits.open(filename+'.fits') as hdul:
	data = hdul[0].data
	hdr = hdul[0].header
	NAXIS1 = hdr['NAXIS1']
	NAXIS2 = hdr['NAXIS2']
	NAXIS3 = hdr['NAXIS3']
	wvl = hdul[1].data[0][0][:,0]
#	wvl = hdul[1].data # for rewitten header
hdulunc = fits.open(uncname+'.fits')
unc = hdulunc[0].data

while True:
	b1 = input("(Re)Analyse? (0 - yes / else - no) ") # keyboard input 1
	if b1=="0":
		wvlnum = int(input("Select a slice of wavelength (0 - {}): ".format(NAXIS3-1))) # keyboard input 2
		print("wavelength = ", wvl[wvlnum])

		## view map 
		figmap = mapview(NAXIS1, NAXIS2, wvl, data, wvlnum, 0)
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

		if not(b0=="2"):
			print("rectangle diagonal ({}, {}) - ({}, {})".format(xb, yb, xt, yt))
			## add patch onto map
			patch = [xb, yb, xt, yt]
			mapview(NAXIS1, NAXIS2, wvl, data, wvlnum, patch)
			## fluxes in box 
			data1 = data[:, yb:yt, xb:xt]
			## view spectrum
			specview(wvl, data1, 0.) #unc_wvl)
			fluxdist(wvl, data1, "all_wvl")
		elif b0=="2":
			print("circle centered at ({0}, {1})".format(x0, y0))
			## add patch onto map
			patch = [x0, y0, radius]
			mapview(NAXIS1, NAXIS2, wvl, data, wvlnum, patch)
			## fluxes in box 
			data1 = [data[:, y0, x0]]
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
					elif (p2!=0 and p2<=r2):
						#print("dx, dy ", dx, dy)
						data1.append(data[:, y0-dy, x0-dx])
						data1.append(data[:, y0-dy, x0+dx])
						data1.append(data[:, y0+dy, x0-dx])
						data1.append(data[:, y0+dy, x0+dx])
			data1 = np.array(data1)
			print("(nb of pix in circle, NAXIS3) = ", data1.shape)
			# view spectrum at selected point
			fig2 = plt.figure("box spectrum")
			ax2 = plt.subplot()
			ax2.errorbar(wvl, np.nansum(data1, axis=0), 0,
				c='k', ls='-', lw=.5, ecolor='r', elinewidth=1., label="IRS data")
			ax2.set_title("Spectrum")
			ax2.set_xlabel('Wavelength (micron)')
			ax2.set_ylabel(r'$F_{\nu} (MJy)$')
#			ax2.set_yscale('log', nonposy='clip')
			ax2.legend(loc='upper left')

		plt.show()
		
		## calculate uncertainty
#		unc_wvl = calunc(wvl, data1, unc1)
		
	else:
		exit(0)
