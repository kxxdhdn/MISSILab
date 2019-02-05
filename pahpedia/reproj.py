import numpy as np
from astropy.io import fits
from reproject import reproject_interp
import hdunit as hd
#from myfunc import deg
import visual as vs


##Search following variables' name to add 'akari'
#chnl
#Nc -> Nc-1
#fits_reproj_filename
#plot_spect_name
#'m83_reproj_0'


fits_image_path = '../data/m83/'
fits_slice_path = '../data/m83/slices/'
fits_conv_path = '../data/m83/convolved/'
fits_reproj_path = '../data/m83/reprojection/'
plot_spect_path = '../plots/m83/spectra/'

##----------choose parametres----------
##channels list
chnl = ['sl2', 'sl1', 'll2', 'll1']
#chnl = ['akari']
Nc = np.size(chnl)
##ref channel for reprojection
fits_ref_filename = 'm83_sl2_10_conv' #demande existence in 'convolved' folder a priori
##-------------------------------------


def rpj(index, wvli, wvl):
	hdr = hd.rhdr(fits_conv_path, fits_ref_filename)[0]
	w = hd.rhdr(fits_conv_path, fits_ref_filename)[15]
	#wvl = []
	flux = []
	k = 0
	##----------boucle_canaux----------
	for i in range(Nc):
		##----------boucle_tranches---------
		for j in range(np.size(wvli[i])):
			if i == Nc: #Nc-1 if add akari
				fits_slice_filename = 'm83_' + chnl[i] + '_' + str(index[i][j])
				with fits.open(fits_slice_path + fits_slice_filename + '.fits') as hdulj:
					dataj, footprint = reproject_interp(hdulj, hdr)#w, shape_out=)
				flux.append(dataj)
			else:
				fits_conv_filename = 'm83_' + chnl[i] + '_' + str(index[i][j]) + '_conv'
				with fits.open(fits_conv_path + fits_conv_filename + '.fits') as hdulj:
					dataj, footprint = reproject_interp(hdulj, hdr)#w, shape_out=)
				flux.append(dataj) #1ere & 2eme dims data
				#print(np.array(flux).shape) ##list dim checker
			k += 1
	data3d = np.array(flux)

	##=======write new fits file=======
	#fits_reproj_filename = 'm83_reproj_0'
	fits_reproj_filename = 'm83_reproj_0'
	comment = 'This FITS file is the homogenized image of m83.'
	hd.cphdr3d(data3d, wvl, 'Wavelength', comment, \
		fits_conv_path, fits_ref_filename, fits_reproj_path, fits_reproj_filename)
	
	##==========plot==========
	spect_title = 'Spectrum'
	for y in range(np.size(data3d[0,:,0])):
		for x in range(np.size(data3d[0,0,:])):
			if all(np.isnan(a)==0 for a in data3d[:,y,x]):
				plot_spect_name = 'm83_spectrum_0' + '_(' + str(x) + ',' + str(y) + ')'
				vs.spect(wvl, data3d, x, y, spect_title, plot_spect_path, plot_spect_name)

	return data3d



def rpj_pp(N_pp, index, wvli):
	hdr = hd.rhdr(fits_conv_path, fits_ref_filename)[0]
	w = hd.rhdr(fits_conv_path, fits_ref_filename)[15]
	flux = []
	k = 0
	##----------boucle_canaux----------
	for i in range(Nc):
		##----------boucle_tranches----------
		for j in range(np.size(wvli[i])):
			if i == Nc: #Nc-1 if add akari
				fits_slice_filename = 'm83_' + chnl[i] + '_' + str(index[i][j])
				with fits.open(fits_slice_path + fits_slice_filename + '.fits') as hdulj:
					dataj, footprint = reproject_interp(hdulj, hdr)#w, shape_out=)
				flux.append(dataj)
			else:
				fits_conv_filename = 'm83_' + chnl[i] + '_' + str(index[i][j]) + '_conv'
				with fits.open(fits_conv_path + fits_conv_filename + '.fits') as hdulj:
					dataj, footprint = reproject_interp(hdulj, hdr)#w, shape_out=)
				flux.append(dataj) #1ere & 2eme dims data
				#print(np.array(flux).shape) ##list dim checker
			k += 1
	data3d = np.array(flux)

	##=======write new fits file=======
	fits_reproj_filename = 'm83_reproj_' + str(N_pp)
	hd.cphdr(data3d, fits_conv_path, fits_ref_filename, fits_reproj_path, fits_reproj_filename)
	
	return data3d
