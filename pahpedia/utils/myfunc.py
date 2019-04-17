#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

My library of math functions

"""

import numpy as np
import subprocess as SP
from .myclass import blockPrint

def fclean(filIN):
	"""
	Clean folder/files
	"""
	SP.call('rm -rf '+filIN, shell=True)
	print("All clean! ")

def hprint(nopr, *content):
	"""
	Hidden print
	"""
	if nopr==True:
		with blockPrint():
			print(*content)
	else:
		print(*content)
"""
def tprint(*content, **kwargs):

	## [TEST] Hidden print

	for key, value in kwargs.items():
		if value==True:
			with blockPrint():
				print([c for c in content])
		else:
			print([c for c in content])
"""
def deg2world(h, m, s, dd, mm, ss):
	"""
	Degree to world coord conversion
	"""
	ra = (h + m/60. + s/3600.) * 360. / 24.
	if dd>0:
		dec = dd + mm/60. + ss/3600.
	else:
		dec = -(-dd + mm/60. + ss/3600.)

	return ra, dec

def gaussian(x, mu, sigma):
	"""
	Normalized Gaussian function given variable x
	"""
	return 1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (x - mu)**2 / (2 * sigma**2))

def gaussian2D(x, y, mux, muy, sigx, sigy, A=1.):
	"""
	2D Gaussian function given iid variables x & y (and amplitude A)
	"""
	return A * np.exp(- (x - mux)**2 / (2 * sigx**2) - (y - muy)**2 / (2 * sigy**2))

def rms(a, ddof=0):
	"""
	Calculate root mean square
	
	--- INPUT ---
	a            an arra or a list
	ddof         Delta Degrees of Freedom
	--- OUTPUT ---
	rms          root mean square of a
	"""
	n = np.size(a) - ddof
	a = np.array(a)
	ms = np.sum(a*a) / n

	return np.sqrt(ms)

def nanrms(a, ddof=0):
	"""
	Calculate root mean square
	
	--- INPUT ---
	a            an arra or a list
	ddof         Delta Degrees of Freedom
	--- OUTPUT ---
	rms          root mean square of a
	"""
	n = np.size(a) - ddof
	a = np.array(a)
	ms = np.nansum(a*a) / n

	return np.sqrt(ms)

def std(a, ddof=0):
	"""
	The same as np.std
	"""
	n = np.size(a) - ddof
	a = np.array(a)
	mu = np.mean(a)
	ms = np.sum((a - mu)**2) / n
	
	return np.sqrt(ms)

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

	print("degree=",deg2world(0, 58, 48, 279, 40, 8))
