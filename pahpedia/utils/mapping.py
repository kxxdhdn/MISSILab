#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

matplotlib applications

"""
from astropy import units as u
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

cmap = mpl.cm.viridis
norm = mpl.colors.Normalize(vmin=0, vmax=1)
COList = ['k', 'b', 'c']
LABEList = ["Spec1", "Spec2", "Spec3"]

def calibview(ref, data, yerr, xerr, figsize=None):
	
	fig, ax = plt.subplots(figsize=figsize)

	ax.plot(np.arange(0,1000,1.), c='grey')
	ax.errorbar(ref, data, yerr=yerr, xerr=xerr, \
		fmt='.', markeredgewidth=.01, capsize=2., c='k', ecolor='r', elinewidth=.5, ls=None)
	ax.set_xlabel("MIPS1 (MJy/sr)")
	ax.set_ylabel("IRS LL1 (MJy/sr)")
	ax.set_xscale('log', nonposx='clip')
	ax.set_yscale('log', nonposy='clip')
	# ax.set_xlim(0, 1e3)
	# ax.set_ylim(0, 1e3)

def specview(wvl, Ldata, wfilt=None, Fnu=None, filt_width=None, figsize=None):

	fig, ax = plt.subplots(figsize=figsize)

	for i, data in enumerate(Ldata):
		ax.errorbar(wvl, data, \
			c=COList[i], ls='-', lw=.5, ecolor='r', elinewidth=1., label=LABEList[i])
			# c=COList[i], ls=None, fmt='.', lw=.5, ecolor='r', elinewidth=1., label=LABEList[i])
		## add inter-calib
		if wfilt!=None:
			photometry = ["IRS/MIPS1", "MIPS1"]
			for j, f in enumerate(Fnu):
				ax.errorbar(wfilt, f, xerr=filt_width,
					fmt='o', capsize=5, c=COList[j+1], ecolor=COList[j+1], elinewidth=.5, label=photometry[j])
	## 
	ax.set_title("Spectra")
	ax.set_xlabel("Wavelength (microns)")
	ax.set_ylabel(r'$F_{\nu} (MJy/sr)$')
	# ax.set_yscale('log', nonposy='clip')
	ax.legend(loc='upper left')
	# fig.savefig('spectrum.png')

def multimview(Larray, w, figshape, abpoints, figsize=None, rot=False, \
	lon_ticksp=1.5 * u.arcmin, lat_ticksp=.5 * u.arcmin, \
	cmap=cmap, norm=norm):
	
	fig, axes = plt.subplots(nrows=figshape[0], ncols=figshape[1], figsize=figsize, \
		subplot_kw=dict(projection=w, sharex=True, sharey=True))

	for i, ax in enumerate(axes.flatten()):
		## define axis
		# ax.set_axis_off()
		## set grid
		ax.coords.grid(color='w', alpha=.5, ls='solid')
		lon = ax.coords[0]
		lat = ax.coords[1]
		# lon.set_major_formatter('hh:mm:ss.s')
		# lat.set_major_formatter('dd:mm:ss.s')
		# lon.set_separator(('h',"'",'"'))
		# lat.set_separator(':')
		lon.set_ticks(spacing=lon_ticksp)
		lat.set_ticks(spacing=lat_ticksp)
		## tick/label position
		if rot==False:
			lat.set_ticks_position('lr')
			lon.set_ticks_position('tb')
			lat.set_ticklabel_position('l')
			lat.set_axislabel_position('l')
			lon.set_ticklabel_position('t')
			lon.set_axislabel_position('t')
			if i==figshape[0]*figshape[1]-1:
				lat.set_ticklabel_position('r')
				lat.set_axislabel_position('r')
				lon.set_ticklabel_position('b')
				lon.set_axislabel_position('b')
		else:
			lat.set_ticks_position('tb')
			lon.set_ticks_position('lr')
			lat.set_ticklabel_position('t')
			lat.set_axislabel_position('t')
			lon.set_ticklabel_position('l')
			lon.set_axislabel_position('l')
			if i==figshape[0]*figshape[1]-1:
				lat.set_ticklabel_position('b')
				lat.set_axislabel_position('b')
				lon.set_ticklabel_position('r')
				lon.set_axislabel_position('r')
		if i==0 or figshape[0]*figshape[1]-1:
			lon.set_axislabel('RA', minpad=.5*3)
			lat.set_axislabel('DEC', minpad=1.*3)
		else:
			## hidden tick label
			lon.set_ticklabel_visible(False)
			lat.set_ticklabel_visible(False)
		## plot
		pcm = ax.imshow(Larray[i], origin='lower', cmap=cmap)
		for pts in abpoints:
			plt.gca()
			plt.plot(pts[0], pts[1], marker='*', c='r')

	fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.85,
                    wspace=0.02, hspace=0.02)

	## add an axes, lower left corner in [0.83, 0.1] measured in figure coordinate with axes width 0.02 and height 0.8
	cb_ax = fig.add_axes([0.88, 0.1, 0.02, 0.8])
	cbar = fig.colorbar(pcm, cax=cb_ax)
	## set the colorbar ticks and tick labels
	# cbar.set_ticks(np.arange(0, 1.1, 0.1))
	# cbar.set_ticklabels(['low', 'medium', 'high'])

"""
------------------------------ MAIN (test) ------------------------------
"""
if __name__ == "__main__":

	import numpy as np
	from rwfits import *

	w = WCSextract('../test_examples/n66_LL1_cube', True)[0]
	multimview(np.array([np.random.random((16,16))]*3), w, (1,3), \
		figsize=(12,8), rot=True, \
		lon_ticksp=.8 * u.arcmin, lat_ticksp=.3 * u.arcmin)

	plt.show()
