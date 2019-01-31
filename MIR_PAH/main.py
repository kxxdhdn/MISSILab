#from os import system
import numpy as np
import math
import hdunit as hd
from myfunc import deg
import decomp as de
import convol as cv
import reproj as rp
import propag as pp


##Search following variables' name to add 'akari'
#'m83_reproj_0'
#fits_reproj_filename
#'wvl.extend(wvli[i])' -> 'wvl = wvli[0]'


fits_image_path = '../data/m83/'
fits_ref0_filename = 'm83_sl2_cube_ref'

w0 = hd.rhdr(fits_image_path, fits_ref0_filename)[15]
NAXIS1, NAXIS2 = hd.rhdr(fits_image_path, fits_ref0_filename)[3:5]

cir = 200


c = input("Do you want to zoom? (YES: enter 'y'/ NO: enter any other key)\n")
if c == 'y':
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
			xmin, ymin = w0.all_world2pix(ramin0, decmin0, 1)
			xmax, ymax = w0.all_world2pix(ramax0, decmax0, 1)
			if xmin > xmax:
				t = xmax
				xmax = xmin
				xmin = t
			if ymin > ymax:
				t = ymax
				ymax = ymin
				ymin = t
			xmin = int(math.floor(xmin))
			xmax = int(math.ceil(xmax))
			ymin = int(math.floor(ymin))
			ymax = int(math.ceil(ymax))
			print('xmin, xmax = ', xmin, xmax)
			print('ymin, ymax = ', ymin, ymax)
			if (0 <= xmin and xmax <= NAXIS1 and 0 <= ymin and ymax <= NAXIS2):
				break
			else:
				print('Invalid number!')
		except ValueError:
			print('Invalid number!')
		except IndexError:
			print('Invalid number!')

##-----test module-----
ramin = ['204', '13', '50']
ramax = ['204', '16', '17']
decmin = ['-29', '-53', '-10']
decmax = ['-29', '-50', '-37']
#ramin = ['204', '14', '52']
#ramax = ['204', '15', '30']
#decmin = ['-29', '-52', '-4']
#decmax = ['-29', '-51', '-34']
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
xmin, ymin = w0.all_world2pix(ramin0, decmin0, 1)
xmax, ymax = w0.all_world2pix(ramax0, decmax0, 1)
if xmin > xmax:
	t = xmax
	xmax = xmin
	xmin = t
if ymin > ymax:
	t = ymax
	ymax = ymin
	ymin = t
xmin = int(math.floor(xmin))
xmax = int(math.ceil(xmax))
ymin = int(math.floor(ymin))
ymax = int(math.ceil(ymax))
print('xmin, xmax = ', xmin, xmax)
print('ymin, ymax = ', ymin, ymax)
##--------------------


##----------without perturbation----------
ii, wvli = de.dec(ramin0, ramax0, decmin0, decmax0)
wvl = []
for i in range(np.size(wvli)):
	wvl.extend(wvli[i])
#wvl = wvli[0]
if isinstance(wvl, list):
	k = 0
	for a in wvl:
		wvl[k] = a
		k += 1
print(np.size(wvl))
#cv.choker(ii, wvli, wvl)
#cv.conv()
#data0 = rp.rpj(ii, wvli, wvl)
#exit(0)
#
###----------'cir' times M-C perturb----------
#flux = []
#for i in range(cir):
#	de.dec_pp(ramin0, ramax0, decmin0, decmax0) ##renew fits files in "slices"
#	cv.conv()
#	data_pp = rp.rpj_pp(i+1, ii, wvli)
#	flux.append(data_pp)
#	print('---------------------------------------------------------------')
#	print(i, '\n')
#	print('---------------------------------------------------------------')


##alternative data reading (data0 & data4)
fits_reproj_path = '../data/m83/reprojection/'
data0, wvl = hd.rdata3d(fits_reproj_path, 'm83_reproj_0')

flux = []
for i in range(cir):
	fits_reproj_filename = 'm83_reproj_' + str(i+1)
	data_pp = hd.rdata(fits_reproj_path, fits_reproj_filename)
	flux.append(data_pp)
##----------------------------------------

data4d = np.array(flux)


##propagation of uncertainty for every pixel
#pp.unc4d(wvl, data0, data4d, cir)
##recadrer reprojected images (in order to pick the pixels which have all wvl data)
#pp.recadre(wvl, data0, data4d, cir)


##Akari-Spitzer
unc = hd.rdata(fits_reproj_path, 'm83_reproj_unc')

data0_a, wvla = hd.rdata3d(fits_reproj_path, 'm83_akari_reproj_0')
unca = hd.rdata(fits_reproj_path, 'm83_akari_reproj_unc')
wvlc = []
data0_c = []
uncc = []
i = np.size(wvla)
while i != 0:
	i -= 1
	wvlc.append(wvla[i])
	data0_c.append(data0_a[i])
	uncc.append(unca[i])
wvlc.extend(wvl)
wvlc = np.array(wvlc)
data0_c.extend(data0)
data0_c = np.array(data0_c)
uncc.extend(unc)
uncc = np.array(uncc)
#exit(0)
pp.spec_comp(wvlc, data0_c, uncc, cir)



