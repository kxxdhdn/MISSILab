#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

matplotlib applications

"""
from astropy import units as u
import numpy as np
from scipy import optimize
import matplotlib as mpl
import matplotlib.pyplot as plt
from myfunc import f_lin, f_lin0
from rwfits import *

cmap = mpl.cm.viridis
norm = mpl.colors.Normalize(vmin=0, vmax=1)
COList = ['k', 'b', 'g', 'm', 'y', 'c']
LABEList = ["Spec1", "Spec2", "Spec3"]

def correview(Rx, Ry, yerr, xerr, source=None, threshold=None, \
	title="Correlation", xlabel="x", ylabel="y", figsize=None, \
	savename='correlation.png'):
	
	fig, ax = plt.subplots(figsize=figsize)
	
	x = np.arange(0,10,.01)

	if source==None:
		source = [input("Input source name: ")]

	for i, s in enumerate(source):

		## linear fit
		fx = Rx[i].ravel()
		fy = Ry[i].ravel()
		valid = ~(np.isnan(fx) | np.isnan(fy))
		popt, pcov = optimize.curve_fit(f_lin, fx[valid], fy[valid])
		# print("popt = ", popt)
		perr = np.sqrt(np.diag(pcov))
		# print("perr = ", perr)

		## draw nstd-sigma intervals
		nstd = threshold
		popt1 = popt + nstd * perr
		popt2 = popt - nstd * perr
		f1 = f_lin(x, *popt1)
		f2 = f_lin(x, *popt2)
		# ax.fill_between(x, f1, f2, facecolor='grey', alpha=.25, \
			# label="{}-sigma intervals".format(threshold))
		ax.errorbar(x, f_lin(x, *popt), \
			c=COList[i], lw=0.5, ecolor='grey', \
			label="{}: Y = {:.2f} * X + {:.2f}".format(s, *popt))
		if s=='M82':
			ax.errorbar(fx, fy, None, None, \
			fmt='.', markeredgewidth=.01, capsize=2., \
			c=COList[i], ecolor=COList[i], elinewidth=.5, ls=None, \
			label=s)
		## scatter
		else:
			ax.errorbar(fx, fy, yerr[i].ravel(), xerr[i].ravel(), \
			fmt='.', markeredgewidth=.01, capsize=2., \
			c=COList[i], ecolor=COList[i], elinewidth=.5, ls=None, \
			label=s)

	ax.set_title(title)
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	ax.set_xscale('log', nonposx='clip')
	ax.set_yscale('log', nonposy='clip')
	ax.set_xlim(.1, 3.)
	ax.set_ylim(.6, 10.)
	ax.legend(loc='upper left')

	fig.savefig(savename)

def calibview(ref, data, yerr, xerr, threshold, \
	xlabel="x", ylabel="y", figsize=None, \
	savename='calib.png'):
	
	fig, ax = plt.subplots(figsize=figsize)

	## y = x
	fx = np.arange(0,1000,1.)
	f1 = fx * threshold
	f2 = fx / threshold
	ax.plot(fx, c='grey')
	ax.fill_between(fx, f1, f2, facecolor='grey', alpha=.25)
	
	## linear fit
	ref = np.array(ref)
	data = np.array(data)
	ref.reshape((-1,))
	data.reshape((-1,))
	valid = ~(np.isnan(ref) | np.isnan(data))
	popt, pcov = optimize.curve_fit(f_lin0, ref[valid], data[valid])
	print("popt = ", popt)
	# perr = np.sqrt(np.diag(pcov))
	# print(perr)
	# ax.errorbar(fx, f_lin0(fx, *popt), yerr=perr[1], c='white', lw=0.5, ecolor='cyan')
	ax.errorbar(fx, f_lin0(fx, *popt), c='b', lw=0.5, ecolor='cyan', \
		label="Y = {:.2f} * X".format(popt[0]))
	factor = 1. / popt[0]
	print("inter-calibration factor (fit) = ", factor)
	
	## (alternative) mean
	r = []
	for i, d in enumerate(data):
		if d!=0 and ref[i]!=0 and np.isnan(d)==0 and np.isnan(ref[i])==0:
			p = ref[i]/d
			if p>.5 and p<2.:
				r.append(p)
	factor2 = np.mean(r)
	print("inter-calibration factor (mean) = ", factor2)

	## scatter
	ax.errorbar(ref, data, yerr=yerr, xerr=xerr, \
		fmt='.', markeredgewidth=.01, capsize=2., c='k', ecolor='r', elinewidth=.5, ls=None)

	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	ax.set_xscale('log', nonposx='clip')
	ax.set_yscale('log', nonposy='clip')
	ax.set_xlim(0, 1e2)
	ax.set_ylim(0, 1e2)
	ax.legend(loc='upperleft')

	fig.savefig(savename)

	return factor

def specview(wvl, Lintensity, Lerr=None, \
	wfilt=None, Fnu=None, filt_width=None, figsize=None, \
	savename='spectrum.png'):

	fig, ax = plt.subplots(figsize=figsize)

	for i, data in enumerate(Lintensity):
		if Lerr==None:
			ax.errorbar(wvl, data, yerr=None, \
				c=COList[i], ls='-', lw=.5, ecolor='r', elinewidth=1., label=LABEList[i])
		else:
			ax.errorbar(wvl, data, yerr=Lerr[i], \
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
	
	fig.savefig(savename)

def imview(Larray, figshape=(1,1), w=gal, add_pts=[], figsize=None, rot=False, \
	lon_ticksp=1.5 * u.arcmin, lat_ticksp=.5 * u.arcmin, \
	cmap=cmap, norm=norm, savename='Image'):
	
	fig, axes = plt.subplots(nrows=figshape[0], ncols=figshape[1], figsize=figsize, \
		subplot_kw=dict(projection=w, sharex=True, sharey=True))
	if figshape==(1,1):
		flaxes = [axes]
	else:
		flaxes = axes.flatten()
	for i, ax in enumerate(flaxes):
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
			if i==0:
				lat.set_ticklabel_position('l')
				lat.set_axislabel_position('l')
				lon.set_ticklabel_position('t')
				lon.set_axislabel_position('t')
			else:
				lat.set_ticklabel_position('r')
				lat.set_axislabel_position('r')
				lon.set_ticklabel_position('b')
				lon.set_axislabel_position('b')
		else:
			lat.set_ticks_position('tb')
			lon.set_ticks_position('lr')
			if i==0:
				lat.set_ticklabel_position('t')
				lat.set_axislabel_position('t')
				lon.set_ticklabel_position('l')
				lon.set_axislabel_position('l')
			else:
				lat.set_ticklabel_position('b')
				lat.set_axislabel_position('b')
				lon.set_ticklabel_position('r')
				lon.set_axislabel_position('r')

		if i==0 or i==figshape[0]*figshape[1]-1:
			lon.set_axislabel('RA', minpad=.5*3)
			lat.set_axislabel('DEC', minpad=1.*3)
		else:
			## hidden tick label
			lon.set_ticklabel_visible(False)
			lat.set_ticklabel_visible(False)
		## plot
		pcm = ax.imshow(Larray[i], origin='lower', cmap=cmap)
		for pts in add_pts:
			plt.gca()
			plt.plot(pts[0], pts[1], marker='*', c='r')

	fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.75,
                    wspace=0.02, hspace=0.02)

	## add an axes, lower left corner in [0.83, 0.1] measured in figure coordinate with axes width 0.02 and height 0.8
	cb_ax = fig.add_axes([0.88, 0.1, 0.02, 0.8])
	cbar = fig.colorbar(pcm, cax=cb_ax)
	## set the colorbar ticks and tick labels
	# cbar.set_ticks(np.arange(0, 1.1, 0.1))
	# cbar.set_ticklabels(['low', 'medium', 'high'])

	fig.savefig(savename)

"""
------------------------------ MAIN (test) ------------------------------
"""
if __name__ == "__main__":

	import numpy as np
	from rwfits import *

	# multimview(np.array([np.random.random((16,16))]*6), \
	# 	figshape=(2,3), figsize=(12,8), rot=True, \
	# 	lon_ticksp=.8 * u.arcmin, lat_ticksp=.3 * u.arcmin)
	imview(np.array([np.random.random((16,16))]*3), \
		figshape=(1,3), figsize=(12,8), rot=True, \
		lon_ticksp=.8 * u.arcmin, lat_ticksp=.3 * u.arcmin)

	plt.show()
