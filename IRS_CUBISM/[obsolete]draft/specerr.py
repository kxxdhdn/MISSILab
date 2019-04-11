from astropy.io import fits
from astropy import units as u
from astropy import wcs
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def spect_res(x, y, wvl, res):
	"""
	plot res - wvl

	--- INPUT ---
	(x, y)     position in pixel coord
	wvl     NAXIS3 elements
	res
	"""
	fig = plt.figure("res distribtion")
	ax = plt.subplot()
	ax.plot(wvl, 0.*wvl, 'k-', linewidth=.5)
	ax.plot(wvl, res, 'b.', markersize=.8)
	ax.set_xlabel("Wavelength (um)")
	ax.set_ylabel("Residual")
	ax.set_ylim(-10, 10)
	
	if np.isnan(x)==0 and np.isnan(y)==0:
		fig.savefig('res/({0}, {1}).png'.format(x, y))
	else:		
		fig.savefig('res/mean.png')

def res_map(x, y, wvl, res):
	"""
	plot residual map

	--- INPUT ---
	x      NAXIS1 elements
	y      NAXIS2 elements
	wvl 1 element
	res
	"""
	fig = plt.figure("res map")
	ax = fig.add_subplot(111, projection='3d')
	ax.scatter(x, y, 0.)
	ax.scatter(x, y, res, s=.5)
	ax.set_xlabel("x")
	ax.set_ylabel("y")
	ax.set_zlabel("Res dist")
	ax.set_zlim(-10,10)

	fig.savefig('res/res_map_{}.png'.format(wvl))

def std_map(x, y, std):
	"""
	plot residual standard deviation map

	--- INPUT ---
	x      NAXIS1 elements
	y      NAXIS2 elements
	std
	"""
	fig = plt.figure("res std dev map")
	ax = fig.add_subplot(111, projection='3d')
	ax.scatter(x, y, 0.)
	ax.scatter(x, y, std, s=.5)
	ax.set_xlabel("x")
	ax.set_ylabel("y")
	ax.set_zlabel("Residual std")
	ax.set_zlim(0,100)

	fig.savefig('res/std_map.png')

def rms(a):
	"""
	calculate root mean square
	
	--- INPUT ---
	a      an array
	"""
	ms = np.nanmean(a*a)

	return np.sqrt(ms)


filename = 'SL1_cube_o'
filename0 = 'SL1_cube_b'
#filename0 = 'sh_SL1_cube'
with fits.open(filename+'.fits') as hdul:
	data = hdul[0].data
with fits.open(filename0+'.fits') as hdul0:
	data0 = hdul0[0].data
	hdr0 = hdul0[0].header
	NAXIS1 = hdr0['NAXIS1']
	NAXIS2 = hdr0['NAXIS2']
	NAXIS3 = hdr0['NAXIS3']
	wvl = hdul0[1].data[0][0][:,0]
res = np.zeros(NAXIS1*NAXIS2*NAXIS3).reshape(NAXIS3, NAXIS2, NAXIS1)
for i in range(NAXIS1):
	for j in range(NAXIS2):
		F_dh = data[:,j,i]
		F_sh = data0[:,j,i]
		for k in range(NAXIS3):
			if np.isnan(F_dh[k])==0 and np.isnan(F_sh[k])==0:
				if F_sh[k]==0.:
					res[k,j,i] = np.nan
				else:
					resk = (F_dh[k] - F_sh[k]) * 2. / (F_dh[k] + F_sh[k])
					res[k,j,i] = resk
			else:
				res[k,j,i] = np.nan

## plot relative residuals distribution at (x, y) or residual (spatial) mean
#spect_res(60, 50, wvl, res[:,50,60])
res_mean = []
for k in range(NAXIS3):
	res_mean.append(np.nanmean(res[k,:,:]))
spect_res(np.nan, np.nan, wvl, res_mean)

## plot residual (spectral) standard deviation map
xs = [] # generate x and y axis
ys = []
k = 15
res_k = []
res_std = []
for j in range(NAXIS2):
	for i in range(NAXIS1):
		xs.append(i)
		ys.append(j)
		res_k.append(res[k,j,i])
		res_std.append(rms(res[:,j,i]))
#res_map(xs, ys, k, res_k)
std_map(xs, ys, res_std)
# statistics
trim = 1.
c = 0
notnan = 0
NPIX = NAXIS1 * NAXIS2
for m in range(NPIX):
	if np.isnan(res_std[m])==0:
		notnan += 1
		if res_std[m]<trim:
			c += 1
print("no NaN percent = ", float(notnan/NPIX))
print("res<{0} percent = {1}".format(trim, float(c/NPIX)))

plt.show()
