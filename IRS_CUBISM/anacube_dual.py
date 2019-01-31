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
def bimapview(x, y, wvl, data, data0, wvlnum, length):
	"""
	view the map and choose 1 point/region to check spectrum & flux distribution

	--- INPUT ---
	x            array to create footprint
	y
	wvl      
	data     
	length       box side length, int type
	--- OUTPUT ---
	(x0, y0)     chosen point
	data1        cut data
	data2        cut data0
	"""

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

	## left click to choose 1 point
	pts = plt.ginput(n=1, show_clicks=True)
	(x0, y0) = np.array(pts[0]).astype(int)
	## fluxes in a square of "length" size centered at the chosen point 
	data1 = data[:, y0-int((length-1)/2):y0+int(length/2)+1, x0-int((length-1)/2):x0+int(length/2)+1]
	data2 = data0[:, y0-int((length-1)/2):y0+int(length/2)+1, x0-int((length-1)/2):x0+int(length/2)+1]
	print("Pixels centered (lowerleft) at ({0}, {1})".format(x0, y0))

	## plot spectrum at chosen square
	fig2 = plt.figure("box spectrum")
	ax2 = plt.subplot()
	ax2.plot(wvl, np.nansum(data1, axis=(1,2)), 'c-', linewidth=.5, label='dh')
	ax2.plot(wvl, np.nansum(data2, axis=(1,2)), 'm-', linewidth=.5, label='sh')
	ax2.set_title("Spectra")
	ax2.set_xlabel('Wavelength (micron)')
	ax2.set_ylabel(r'$F_{\nu} (MJy)$')
#	ax2.set_yscale('log', nonposy='clip')
	ax2.legend(loc='best')
#	fig2.savefig('spectra2.png')

	"""
	box statistics
	"""
	## plot flux distribution for pixels in the box
	fig3 = plt.figure("box flux distribution")
	ax3 = plt.subplot()
	# ----- target -----
	num, bins, patches = ax3.hist(data1[wvlnum,:,:].reshape(-1), bins=int(length*length/4),
		alpha=.3, lw=.5, edgecolor='w', color='c')
	bins_c = []
	for i in range(np.size(bins)-1):
		bins_c.append((bins[i]+bins[i+1])/2.)
	bins_c = np.array(bins_c)
	# ----- ref -----
	num0, bins0, patches0 = ax3.hist(data2[wvlnum,:,:].reshape(-1), bins=int(length*length/4),
		alpha=.3, lw=.5, edgecolor='w', color='m')
	bins0_c = []
	for i in range(np.size(bins0)-1):
		bins0_c.append((bins0[i]+bins0[i+1])/2.)
	bins0_c = np.array(bins0_c)

	## Fit the data using a Gaussian
	g_init = models.Gaussian1D(amplitude=1., mean=0, stddev=1.)
	fit_g = fitting.LevMarLSQFitter()
	# ----- target -----
	g = fit_g(g_init, bins_c, num)
	numsum = 0
	for j in range(np.size(num)):
		if (bins_c[j]>g.mean-g.stddev and bins_c[j]<g.mean+g.stddev):
			numsum += num[j]
	sig1 = float(numsum / (length*length))
	print("box flux mean: ", g.mean*1.)
	print("1-sigma percentage: {:.0%}".format(sig1))
	ax3.plot(bins_c, g(bins_c), c='c')
	ax3.axvline(g.mean-g.stddev, c='c', linestyle=':')
	ax3.axvline(g.mean+g.stddev, c='c', linestyle=':')
	# ----- ref -----
	g0 = fit_g(g_init, bins0_c, num0)
	num0sum = 0
	for j in range(np.size(num0)):
		if (bins0_c[j]>g0.mean-g0.stddev and bins0_c[j]<g0.mean+g0.stddev):
			num0sum += num0[j]
	sig1ref = float(num0sum / np.nansum(num))
	print("ref box flux mean: ", g0.mean*1.)
	print("ref 1-sigma percentage: {:.0%}".format(sig1ref))
	ax3.plot(bins0_c, g0(bins0_c), c='m')
	ax3.axvline(g0.mean-g0.stddev, c='m', linestyle=':')
	ax3.axvline(g0.mean+g0.stddev, c='m', linestyle=':')

	ax3.set_title("Flux distribution")
	ax3.set_xlabel(r'$F_{\nu} (MJy)$')
	ax3.set_ylabel("Number of pixels")

	plt.show()

	return x0, y0, data1, data2

def cal_res(wvl, data1, data2, length):
	"""
	calculate residuals for chosen pixels with each wavelength

	--- INPUT ---
	wvl       
	data1     
	data2     
	length       int type
	--- OUTPUT ---
	res          3D array
	"""
	NX = length
	NY = length
	NZ = np.size(wvl)
	res = np.zeros(NX*NY*NZ).reshape(NZ, NY, NX) # create a 3D array
	for i in range(NX):
		for j in range(NY):
			F_dh = data1[:,j,i]
			F_sh = data2[:,j,i]
			for k in range(NZ):
				if np.isnan(F_dh[k])==0 and np.isnan(F_sh[k])==0:
					if F_sh[k]==0.:
						res[k,j,i] = np.nan
					else:
						resk = (F_dh[k] - F_sh[k]) * 2. / (F_dh[k] + F_sh[k])
						res[k,j,i] = resk
				else:
					res[k,j,i] = np.nan

	return res

def spec_res(x0, y0, wvl, res):
	"""
	plot res - wvl & res_std. - wvl

	--- INPUT ---
	(x0, y0)     position in pixel coord
	wvl          NAXIS3 elements
	res          3D array
	"""
	fig = plt.figure("res & res_std distribution", figsize=(14,5))
	## res - wvl
	res_wvl = np.nanmean(res, axis=(1,2))
	res_wvl_mean = np.nanmean(res_wvl)
	print("res_wvl_mean = ", res_wvl_mean)
	ax1 = plt.subplot(121)
	ax1.plot(wvl, 0.*wvl, 'k-', linewidth=.5)
#	ax1.plot(wvl, res_wvl_mean, 'k-', linewidth=.5)
	ax1.plot(wvl, res_wvl, 'b.', markersize=.8)
	ax1.set_xlabel("Wavelength (um)")
	ax1.set_ylabel("Residual")
	ax1.set_ylim(-3, 3)
	## res_std - wvl
	std_wvl = []
	for k in range(np.size(wvl)):
		std_wvl.append(rms(res[k,:,:]))
	std_wvl_mean = np.nanmean(std_wvl)
	print("std_wvl_mean = ", std_wvl_mean)
	ax2 = plt.subplot(122)
	ax2.plot(wvl, 0.*wvl, 'k-', linewidth=.5)
#	ax2.plot(wvl, std_wvl_mean, 'k-', linewidth=.5)
	ax2.plot(wvl, std_wvl, 'b.', markersize=.8)
	ax2.set_xlabel("Wavelength (um)")
	ax2.set_ylabel("Residual std")
	ax2.set_ylim(-3, 3)

#	fig.savefig('res/({0}, {1}).png'.format(x0, y0))
	plt.show()

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
filename = 'dh_SL2_cube'
filename0 = 'sh_SL2_cube' 
#filename0 = filename # set 2 filenames the same if no need to compare cubes

with fits.open(filename+'.fits') as hdul:
	data = hdul[0].data
	hdr = hdul[0].header
	NAXIS1 = hdr['NAXIS1']
	NAXIS2 = hdr['NAXIS2']
	NAXIS3 = hdr['NAXIS3']
	wvl = hdul[1].data[0][0][:,0]
#	wvl = hdul[1].data # for rewitten header
with fits.open(filename0+'.fits') as hdul0:
	data0 = hdul0[0].data

while True:
	flag = input("(Re)Analyse? (0 - yes / else - no) ") # keyboard input 1
	if (flag=="0"):
		wvlnum = int(input("Choose a slice of wavelength: ")) # keyboard input 2
		print("wavelength = ", wvl[wvlnum])
		length = int(input("Choose the square length: ")) # keyboard input 3
		x = np.arange(hdr['NAXIS1'])
		y = np.arange(hdr['NAXIS2'])
		x0, y0, data1, data2 = bimapview(x, y, wvl, data, data0, wvlnum, length)
		"""
		compare two cubes [comment the following if not compare]
		"""
#		res = cal_res(wvl, data1, data2, length)
#		print("res shape = ", res.shape)
#		spec_res(x0, y0, wvl, res)
	else:
		exit(0)
