import numpy as np

def rms(x):
	n = np.size(x) - 1
	x = np.array(x)
	mu = np.mean(x)
	ms = np.sum((x - mu)**2)
	y = np.sqrt(ms / n)
	return y



def gaussian(x, mu, sigma):
	return 1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (x - mu)**2 / (2 * sigma**2))



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



def deg(dd, mm, ss):
	d = dd + mm/60. + ss/3600.
	return d