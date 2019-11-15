#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

My library of useful functions

"""

import math
import numpy as np
import scipy.interpolate as interpolate
import subprocess as SP

def fclean(filIN, *alert):
	'''
	Clean folder/files
	'''
	SP.call('rm -rf '+filIN, shell=True)
	for text in alert:
		print(text)

def pix2sr(X, CDELT):
	'''
	X pixel = Y sr
	'''
	PFOV = 3600 * abs(CDELT)
	
	return X * (PFOV * 2. * math.pi / (360. * 3600.))**2.

def sr2arcsec2(X):
	"""
	X sr = Y arcsec^2
	"""
	return X * (360. * 3600. / (2. * math.pi))**2.

def rad2arcsec(X):
	'''
	X rad = Y arcsec
	'''
	return X * 360. * 3600. / (2. * math.pi)

def celest2deg(h, m, s, dd, mm, ss):
	'''
	Degree to world coord conversion
	'''
	ra = (h + m/60. + s/3600.) * 360./24.
	if dd>0:
		dec = dd + mm/60. + ss/3600.
	else:
		dec = -(-dd + mm/60. + ss/3600.)

	return ra, dec

def deg2celest(ra, dec):
	'''
	World coord to degree conversion
	'''
	h = math.floor(ra * 24./360.)
	m = math.floor((ra - h) * 60.)
	s = (ra - h - m/60.) * 3600.
	if dec>0:
		dd = math.floor(dec)
		mm = math.floor((dec - dd) * 60.)
		ss = (dec - dd - mm/60.) * 3600.
	else:
		dd = math.ceil(dec) # negtive
		mm = -math.ceil((dec - dd)) # positive
		ss = -(dec - dd + mm/60.) * 3600. # positive

	print('{.0f}h{.0f}m{04.2f}s, {.0f}d{.0f}m{04.2f}s'.format(h,m,s,dd,mm,ss))

	return h, m, s, dd, mm, ss

def f_lin(x, A, B):
	'''
	Y = A * x + B
	'''
	return A * x + B

def f_lin1(x, B):
	'''
	Y = x + B
	'''
	return x + B

def f_lin0(x, A):
	'''
	Y = A * x
	'''
	return A * x

def gaussian(x, mu, sigma):
	'''
	Normalized Gaussian function given variable x
	'''
	return 1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (x - mu)**2 / (2 * sigma**2))

def gaussian2D(x, y, mux, muy, sigx, sigy, A=1.):
	'''
	2D Gaussian function given iid variables x & y (and amplitude A)
	'''
	return A * np.exp(- (x - mux)**2 / (2 * sigx**2) - (y - muy)**2 / (2 * sigy**2))

def rms(a, ddof=0):
	'''
	Calculate root mean square
	
	--- INPUT ---
	a            an array or a list
	ddof         Delta Degrees of Freedom
	--- OUTPUT ---
	rms          root mean square of a
	'''
	n = np.size(a) - ddof
	a = np.array(a)
	ms = np.sum(a*a) / n

	return np.sqrt(ms)

def nanrms(a, ddof=0):
	'''
	Calculate root mean square (nan = 0)
	
	--- INPUT ---
	a            an array or a list
	ddof         Delta Degrees of Freedom
	--- OUTPUT ---
	rms          root mean square of a
	'''
	n = np.size(a) - ddof
	a = np.array(a)
	ms = np.nansum(a*a) / n

	return np.sqrt(ms)

def std(a, ddof=0):
	'''
	The same as np.std
	'''
	n = np.size(a) - ddof
	a = np.array(a)
	mu = np.mean(a)
	ms = np.sum((a - mu)**2) / n
	
	return np.sqrt(ms)

def closest(a, val):
	'''
	Return the index i corresponding to the closest a[i] to val
	'''
	a = list(a)
	
	return a.index(min(a, key=lambda x:abs(x-val)))

def bsplinterpol(x, y, x0):
	'''
	Monte-Carlo propagated error calculator
	
	--- INPUT ---
	x           in base x
	y           in data y
	x0          out base x
	--- OUTPUT ---
	bspl(x0)    B-spline interpol out data
	'''
	mask = []
	for i, yi in enumerate(y):
		if np.isnan(yi)==1 or yi==0:
			mask.append(i)
	if len(x)-len(mask)>4: # number of knots (avoid all NaNs col)
		x = np.delete(x, mask)
		y = np.delete(y, mask)

	t, c, k = interpolate.splrep(x, y, s=0, k=4) # s - smooth
	# print('''\
	# t: {}
	# c: {}
	# k: {}
	# '''.format(t, c, k))
	bspl = interpolate.BSpline(t, c, k, extrapolate=False)
	
	return bspl(x0)

def MCerror(arr, axis=0):
	'''
	Monte-Carlo propagated error calculator
	
	--- INPUT ---
	arr         array
	axis        axis or axes along which error is calculated
	--- OUTPUT ---
	err         one dim less than arr
	'''
	## Detect dimension
	##------------------
	axsh = arr.shape
	NAXIS = np.size(axsh)
	dim = NAXIS - 1
	if axis==0:
		if dim==2:
			err = np.full((axsh[1],axsh[2]), np.nan)
			for i in range(axsh[2]):
				for j in range(axsh[1]):
					err[j,i] = np.nanstd(arr[:,j:i], ddof=1)
		elif dim==3:
			err = np.full((axsh[1],axsh[2],axsh[3]), np.nan)
			for i in range(axsh[3]):
				for j in range(axsh[2]):
					for k in range(axsh[1]):
						err[k,j,i] = np.nanstd(arr[:,k,j,i], ddof=1)
		else:
			print('ERROR: dimension not supported! ')
	else:
		print('ERROR: array shape not supported! ')

	return err

"""
------------------------------ MAIN (test) ------------------------------
"""
if __name__ == "__main__":

	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D

	#x=np.random.random(100)
	x=np.arange(-5.,5.,.1)
	#print("x=", x)
	
	std0=np.std(x,ddof=0)
	print("std by np =",std0)
	std=std(x,ddof=0)
	print("std=",std)
	
	rms=rms(x,ddof=0)
	print("rms=",rms)
	
	g=gaussian(x, 0., 1.)
	#print(g)
	plt.plot(x,g)
	g2=gaussian2D(x, x, 0., 0., 1., 1.)
	fig = plt.figure("res std dev map")
	ax = fig.add_subplot(111, projection='3d')
	ax.scatter(x, x, g2, s=.5)
	#plt.show()

	print("degree=",celest2deg(0, 58, 48, -179, 40, 8))
