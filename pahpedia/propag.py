import numpy as np
import hdunit as hd
from myfunc import rms
import visual as vs


##Search following variables' name to add 'akari'
#plot_spec_name
#fits_err_filename
#'m83_reproj_0'
#fits_reproj_filename
#vs.spect_err(data -> data0)


plot_per_path = '../plots/m83/perturbation/'
plot_spec_path = '../plots/m83/spectra/'

fits_reproj_path = '../data/m83/reprojection/'
fits_conv_path = '../data/m83/convolved/'
fits_ref_filename = 'm83_sl2_10_conv'
fits_slice_path = '../data/m83/slices/'


def unc4d(wvl, data0, data4d, N_pp):
	w = hd.rhdr(fits_conv_path, fits_ref_filename)[15]
	data = hd.rdata(fits_reproj_path, 'm83_reproj_brute')
	NAXISx = np.size(data4d[0,0,0,:])
	NAXISy = np.size(data4d[0,0,:,0])
	NAXISz = np.size(data4d[0,:,0,0])
	im_err = np.full((NAXISz, NAXISy, NAXISx), np.nan)
	for y in range(NAXISy):
		for x in range(NAXISx):
#			##----------plot distribution----------
			if all(np.isnan(a)==0 for a in data4d[0,:,y,x]):
				sig_per = []
				for i in range(N_pp):
					per = []
					for a3 in data4d[i,:,y,x]: #
						per.append(a3)
					sig_per.append(rms(per))
				plot_per_name = 'per_spec_(' + str(x) + ',' + str(y) + ')'
				vs.dist(sig_per, plot_per_path, plot_per_name)
			##----------plot spectra----------
			for z in range(NAXISz):
				if all(np.isnan(a)==0 for a in data4d[:,z,y,x]):
					per = []
					for a4 in data4d[:,z,y,x]: #return 'yerr' with a size of z (nb of wvl)
						per.append(a4)
					im_err[z,y,x] = rms(per)
			if all(np.isnan(a)==0 for a in data0[:,y,x]):
				spect_title = 'Spectrum of M83'
				plot_spec_name = 'm83_spitzer_spectrum_(' + str(x) + ',' + str(y) + ')'
				vs.spect_err(wvl, data[250:], data0, x, y, im_err[:,y,x], spect_title, plot_spec_path, plot_spec_name)
	##----------save yerr----------
#				##calculate S/N
				noise = np.mean(im_err[:,y,x])
				signal = np.mean(data0[:,y,x])
				sn = signal/noise
				print(sn)
	fits_err_filename = 'm83_reproj_unc'
	comment = 'This is the uncertainty file of the homogenized image of M83.'
	hd.cphdr3d(im_err,  wvl, 'Wavelength', comment, \
		fits_conv_path, fits_ref_filename, fits_reproj_path, fits_err_filename)



def spec_comp(wvl, data0, unc, N_pp):
	w = hd.rhdr(fits_conv_path, fits_ref_filename)[15]
	##Reprojection without perturbation
	data = hd.rdata(fits_reproj_path, 'm83_reproj_brute')
	NAXISx = np.size(data0[0,0,:])
	NAXISy = np.size(data0[0,:,0])
	NAXISz = np.size(data0[:,0,0])
	for y in range(NAXISy):
		for x in range(NAXISx):
			if all(np.isnan(a)==0 for a in data0[:,y,x]):
				print('(', x, ',', y, ')')
			if any(np.isnan(a)==0 for a in data0[:,y,x]):
				spect_title = 'Spectrum of M83'
				plot_spec_name = 'm83_spectrum_(' + str(x) + ',' + str(y) + ')'
				vs.spect_err_c(wvl, data, data0, x, y, unc[:,y,x], \
					spect_title, plot_spec_path, plot_spec_name)



def recadre(wvl, data0, data4d, N_pp):
	w = hd.rhdr(fits_conv_path, fits_ref_filename)[15]
	NAXISx = np.size(data4d[0,0,0,:])
	NAXISy = np.size(data4d[0,0,:,0])
	NAXISz = np.size(data4d[0,:,0,0])
	##----------rewrite reproj images----------
	##-----reproj_0-----
	newdata = np.full((NAXISz, NAXISy, NAXISx), np.nan)
	for y in range(NAXISy):
		for x in range(NAXISx):
			if all(np.isnan(a)==0 for a in data0[:,y,x]):
				for z in range(NAXISz):
					newdata[z,y,x] = data0[z,y,x]
	#comment = 'This FITS file is the homogenized image of M83.'
	#hd.cphdr3d(newdata, wvl, 'Wavelength', comment, \
	#	fits_conv_path, fits_ref_filename, fits_reproj_path, 'm83_reproj_0')
	##-----reproj_i-----
	for i in range(N_pp):
		newdata = np.full((NAXISz, NAXISy, NAXISx), np.nan)
		for y in range(NAXISy):
			for x in range(NAXISx):
				if all(np.isnan(a)==0 for a in data4d[i,:,y,x]): #Attention! instead of i+1
					for z in range(NAXISz):
						newdata[z,y,x] = data4d[i,z,y,x]
		#comment = 'This FITS file is perturbed image NO.' + str(i+1) + 'of M83.'
		#fits_reproj_filename = 'm83_reproj_' + str(i+1)
		#hd.cphdr3d(newdata, wvl, 'Wavelength', comment, \
		#	fits_conv_path, fits_ref_filename, fits_reproj_path, fits_reproj_filename)

