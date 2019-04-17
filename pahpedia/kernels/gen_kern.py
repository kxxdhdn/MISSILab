#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Kernel generation

"""

import numpy as np

def kern(lam):
	sim_par_wave = [0, 13.25, 40.]
	sim_par_fwhm = [2.8, 3.26, 10.1]
	sim_per_wave = [0, 15.5, 40.]
	sim_per_fwhm = [3.8, 3.8, 10.1]
	##fwhm (arcsec)
	fwhm_par = np.interp(lam, sim_par_wave, sim_par_fwhm)
	fwhm_per = np.interp(lam, sim_per_wave, sim_per_fwhm)
	#fwhm_lam = np.sqrt(fwhm_par * fwhm_per)
	##sigma (arcsec)
	sim_par = fwhm_par / (2. * np.sqrt(2.*np.log(2.)))
	sim_per = fwhm_per / (2. * np.sqrt(2.*np.log(2.)))
	sim_lam = np.sqrt(sim_par * sim_per)
	return sim_lam, sim_par, sim_per

"""
------------------------------ MAIN (test) ------------------------------
"""
if __name__ == "__main__":

	import matplotlib.pyplot as plt
	
	path = '../test_examples/'
	wvl = np.arange(5., 40., .25)
	sigma_lam, sigma_par, sigma_per = kern(wvl)
	fig = plt.figure()
	plt.plot(wvl, sigma_par, color='g', label='Parallel')
	plt.plot(wvl, sigma_per, color='c', label='Perpendicular')
	plt.plot(wvl, sigma_lam, color='k', label='Geometric mean')
	plt.title('PSF profile')
	plt.xlabel('Wavelength, $\lambda$ [$\mu$m]')
	plt.ylabel('PSF [arcsec]')
	plt.legend(loc='upper left')
	
	plt.savefig(path + 'gen_kern.jpeg')
	#plt.show()