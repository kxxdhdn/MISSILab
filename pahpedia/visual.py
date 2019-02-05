import numpy as np
import math
import astropy.units as u
import matplotlib as mpl
import matplotlib.pyplot as plt
from myfunc import deg, gaussian


def pj_wcs(data, w, arcs, im_title):
	plt.figure(figsize=(3, 6))
	cmap = mpl.cm.viridis
	norm = mpl.colors.Normalize(vmin=0, vmax=1)

	ax = plt.subplot(projection=w)
	pcm = ax.imshow(data, origin='lower', cmap=cmap)
	plt.colorbar(pcm, cmap=cmap)
	coord = ax.coords
	coord.grid(color='y', alpha=0.5)
	lon = coord[0]
	lon.set_axislabel('Right Ascension', minpad=1)
	lon.set_major_formatter('dd:mm:ss.s')
	lon.set_separator(('°', "'", '"'))
	lon.set_ticks(spacing=arcs*u.arcsec)
	lat = coord[1]
	lat.set_axislabel('Declination', minpad=1)
	lat.set_major_formatter('dd:mm:ss.s')
	lat.set_separator(('°', "'", '"'))
	lat.set_ticks(spacing=arcs*u.arcsec)
	
	nx = np.size(data[0,:])
	ny = np.size(data[:,0])
	print('Image size in pixel: ', nx, '*', ny)
	##choose point?
	c1 = input("Do you want to choose a point? (YES: enter 'y'/ NO: enter any other key)\n")
	if c1 == 'y':
		while True:
			ra = input("RA = (dd mm ss) ").split()
			dec = input("DEC = (dd mm ss) ").split()
			try:
				for i in range(3):
					ra[i] = float(ra[i])
					dec[i] = float(dec[i])
				ra0 = deg(ra[0], ra[1], ra[2])
				dec0 = deg(dec[0], dec[1], dec[2])
				## WCS - pixel coords convert
				x0, y0 = w.all_world2pix(ra0, dec0, 1)
				#ax.plot(x0, y0, '+', color='r')
				x0 = int(math.floor(x0))
				y0 = int(math.floor(y0))
				print('x0, y0 = ', x0, y0)
				if (0 <= x0 <= nx and 0 <= y0 <= ny):
					break
				else:
					print('Invalid number!')
			except ValueError:
				print('Invalid number!')
			except IndexError:
				print('Invalid number!')
		ax.plot(x0, y0, '+', color='r')

	##zoom?
	c2 = input("Do you want to zoom? (YES: enter 'y'/ NO: enter any other key)\n")
	if c2 == 'y':
		while True:
			ramin = input("RA_min = (dd mm ss) ").split()
			ramax = input("RA_max = (dd mm ss) ").split()
			decmin = input("DEC_min = (dd mm ss) ").split()
			decmax = input("DEC_max = (dd mm ss) ").split()
			try:
				for i in range(3):
					ramin[i] = float(ramin[i])
					ramax[i] = float(ramax[i])
					decmin[i] = float(decmin[i])
					decmax[i] = float(decmax[i])
				ramin0 = deg(ramin[0], ramin[1], ramin[2])
				ramax0 = deg(ramax[0], ramax[1], ramax[2])
				decmin0 = deg(decmin[0], decmin[1], decmin[2])
				decmax0 = deg(decmax[0], decmax[1], decmax[2])
				## WCS - pixel coords convert
				xmin, ymin = w.all_world2pix(ramin0, decmin0, 1)
				xmax, ymax = w.all_world2pix(ramax0, decmax0, 1)
				if xmin > xmax:
					t = xmax
					xmax = xmin
					xmin = t
				if ymin > ymax:
					t = ymax
					ymax = ymin
					ymin = t
				xmin = math.floor(xmin)
				xmax = math.ceil(xmax)
				ymin = math.floor(ymin)
				ymax = math.ceil(ymax)
				if (0 <= xmin and xmax <= nx and 0 <= ymin and ymax <= ny):
					break
				else:
					print('Invalid number!')
			except ValueError:
				print('Invalid number!')
			except IndexError:
				print('Invalid number!')
		ax.set_xlim(xmin, xmax)
		ax.set_ylim(ymin, ymax)
	
	ax.set_title(im_title)
	#ax.autoscale(False)

	##save plot?
	c3 = input("Do you want to save the plot? (YES: enter 'y'/ NO: enter any other key)\n")
	if c3 == 'y':
		plot_save_path = input("Path = ")
		plot_save_name = input("Name = ")
		plt.savefig(plot_save_path+plot_save_name+'.jpeg')

	plt.show()



def spect(wvl, data3d, x0, y0, im_title, plot_path, plot_name):
	plt.figure()

	ax = plt.subplot()
	ax.plot(wvl, data3d[:,y0,x0], color='r')
	ax.set_title(im_title)
	#ax.set_xlim(7,8)
	#ax.set_xlim(14,15)
	#ax.set_xlim(20,21)

	plt.savefig(plot_path+plot_name+'.jpeg')
	##please don't plt.show()



def spect_err(wvl, data, data3d, x0, y0, yerr, im_title, plot_path, plot_name):
	plt.figure()

	ax = plt.subplot()
	ax.plot(wvl, data[:,y0,x0], color='y')
	ax.errorbar(wvl, data3d[:,y0,x0], yerr, c='k', ls='-', lw=0.5, ecolor='r', elinewidth=1.5)
	ax.set_title(im_title)
	ax.set_yscale('log', nonposy='clip')
	#ax.set_xlim(6,8.5)

	plt.savefig(plot_path+plot_name+'.jpeg')
	##please don't plt.show()



def spect_err_c(wvl, data, data3d, x0, y0, yerr, im_title, plot_path, plot_name):
	plt.figure(figsize=(8,5))

	ax = plt.subplot()
	ax.plot(wvl[:250], data[:250,y0,x0], color='y')
	ax.plot(wvl[252:], data[252:,y0,x0], color='y')
	ax.errorbar(wvl[:250], data3d[:250,y0,x0], yerr[:250], c='k', ls='-', lw=0.5, ecolor='r', elinewidth=1.5)
	ax.errorbar(wvl[252:], data3d[252:,y0,x0], yerr[252:], c='k', ls='-', lw=0.5, ecolor='r', elinewidth=1.5)
	ax.set_title(im_title)
	ax.set_xlabel('Wavelength, $\lambda$ [$\mu$m]')
	ax.set_ylabel('Flux, $F_{\mu}$ [MJy/sr]')
	ax.set_xscale('log', nonposy='clip')
	ax.set_yscale('log', nonposy='clip')
	#ax.set_xlim(6,8.5)

	plt.savefig(plot_path+plot_name+'.jpeg')
	##please don't plt.show()



def dist(flux, plot_path, plot_name):
	plt.figure()

	count, bins, ignored = plt.hist(flux, bins='auto', \
		alpha=0.5, lw=0.5, edgecolor='k', color='c')#, normed=True)
	plt.title('histogram')
	plt.xlabel('$\sigma_{\mu}$ [MJy/sr]')
	plt.ylabel('Number of pixels')
	#plt.plot(bins, gaussian(bins,  mu, sigma), lw=2, color='r')

	plt.savefig(plot_path+plot_name+'.jpeg')
	##please don't plt.show()
