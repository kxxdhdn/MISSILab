from os import system
import numpy as np
import math
import hdunit as hd
import visual as vs
import decomp as dc
import reproj as rp
import convol as cv
from myfunc import deg, rms
#import propag as pp

chnl = ['sl2', 'sl1', 'll2', 'll1']

fits_image_path = "../data/m83/"
fits_image_filename = 'm83_sl2_cube'
fits_ref0_filename = 'm83_sl2_cube_ref'
#fits_image_filename = 'm83_' + chnl[i] + '_cube'
#fits_unc_filename = 'm83_' + chnl[i] + '_cube_unc'
#fits_ref0_filename = 'm83_' + chnl[i] + '_cube_ref'

fits_slice_path = '../data/m83/slices/'
#fits_slice_filename = 'm83_' + chnl[i] + '_' + str(j)

fits_conv_path = '../data/m83/convolved/'
fits_ref_filename = 'm83_sl2_0_conv'

fits_reproj_path = '../data/m83/reprojection/'
#fits_reproj_filename = 'm83_reproj_' + str(N_circle)

plot_save_path = '../plots/m83/'
plot_save_name = 'test'


#NAXIS1, NAXIS2, NAXIS3, CRPIX1, CRPIX2 = hd.rhdr3d(fits_image_path, fits_image_filename)[3:8]
#ax = np.arange(NAXIS1) #1--ax
#ay = np.arange(NAXIS2) #2--ay

##WCS reprojection
#data = hd.rdata3d(fits_image_path, fits_image_filename)[0]
#data = hd.rdata(fits_image_path, 'm83_ll2_cube')[10]
#w = hd.rhdr(fits_image_path, 'm83_ll2_cube_ref')[15]
#data = hd.rdata(fits_slice_path, 'm83_ll2_10')
#w = hd.rhdr(fits_slice_path, 'm83_ll2_10')[15]
#vs.pj_wcs(data, w, 30., 'akari100')

#data = hd.rdata3d(fits_reproj_path, 'm83_reproj_brute')[0]
#print(data.shape)
#w = hd.rhdr(fits_conv_path, 'm83_ll2_10_conv')[15]
#vs.pj_wcs(data[190,:,:], w, 30., 'M83 reprojection')

data = hd.rdata(fits_image_path, 'm83_akari_cube')
print(data.shape)
w = hd.rhdr(fits_image_path, 'm83_akari_cube_ref')[15]
vs.pj_wcs(data[23,:,:], w, 5., 'AKARI FOV (M83)')

##plot
#vs.pj_wcs(data[20,:,:], w, 30., 'test_wcs')#, 117, 148, 117, 129)

#hd.cphdr(data[0,:,:], fits_image_path, fits_image_filename, plot_save_path, plot_save_name)
#hd.cphdr(data, fits_image_path, fits_image_filename, plot_save_path, 'test')


##----------rewrite brute.fits file----------
#fits_conv_path = '../data/m83/convolved/'
#fits_ref_filename = 'm83_sl2_10_conv'
#
#data0, wvl0 = hd.rdata3d(fits_image_path, 'm83_reproj_brute')
#dataa, wvla = hd.rdata3d(fits_reproj_path, 'm83_akari_reproj_0')
#
#wvl=[]
#data=[]
#i = 250
#while i != 0:
#	i -= 1
#	wvl.append(wvla[i])
#	data.append(dataa[i])
#for i in range(344):
#	wvl.append(wvl0[i])
#	data.append(data0[i])
#data = np.array(data)
#print(data.shape)
#comment = 'This FITS file is the non-homogenized image of m83.'
#hd.cphdr3d(data, wvl, 'Wavelength', comment, \
#		fits_conv_path, fits_ref_filename, '../data/m83/', 'm83_reproj_brute')
##------------------------------


##----------rebuilt akari unc----------
##-----choose region-----
#w0 = hd.rhdr(fits_image_path, 'm83_akari_cube_ref')[15]
#NAXIS1, NAXIS2 = hd.rhdr(fits_image_path, 'm83_akari_cube_ref')[3:5]
#w = hd.rhdr(fits_slice_path, 'm83_ll2_10')[15]
#
#ramin0 = 204.2617
#ramax0 = 204.2637
#decmin0 = -29.8672
#decmax0 = -29.8660
### WCS - pixel coords convert
#xmin, ymin = w0.all_world2pix(ramin0, decmin0, 1)
#xmax, ymax = w0.all_world2pix(ramax0, decmax0, 1)
#if xmin > xmax:
#	t = xmax
#	xmax = xmin
#	xmin = t
#if ymin > ymax:
#	t = ymax
#	ymax = ymin
#	ymin = t
#xmin = int(math.floor(xmin))
#xmax = int(math.ceil(xmax))
#ymin = int(math.floor(ymin))
#ymax = int(math.ceil(ymax))
#print('xmin, xmax = ', xmin, xmax)
#print('ymin, ymax = ', ymin, ymax)
###--------------------
#data = hd.rdata(fits_image_path, 'm83_akari_cube')[21:271]
#wvl = hd.rdata3d(fits_reproj_path, 'm83_akari_reproj_0')[1]
#data = data[:, ymin:ymax+1, xmin:xmax+1]
#vs.pj_wcs(data[0,:,:], w, 30., 'akari100')
#unc = np.full((np.size(wvl), NAXIS2, NAXIS1), np.nan)
#for i in range(np.size(wvl)):
#	datai = []
#	for a in data[i,:,:].reshape((ymax-ymin+1)*(xmax-xmin+1)):
#		if np.isnan(a)==0:
#			datai.append(a)
#	datai = np.array(datai)
#	unc[i,:,:] = rms(datai)
#	#print(rms(datai))
##print(unc[101,:,:])
#
##vs.pj_wcs(unc[0], w, 30., 'akari100')
#hd.cphdr(unc, fits_image_path, 'm83_akari_cube_ref', fits_image_path, 'm83_akari_cube_unc')
#hd.cphdr(data, fits_image_path, 'm83_akari_cube_ref', fits_image_path, 'test')
##------------------------------

