import numpy as np
import math
import hdunit as hd


##Search following variables' name to add 'akari'
#chnl
#Nc -> Nc-1


fits_image_path = '../data/m83/'
fits_slice_path = '../data/m83/slices/'
chnl = ['sl2', 'sl1', 'll2', 'll1']
#chnl = ['akari']
Nc = np.size(chnl)


def dec(ramin, ramax, decmin, decmax): #dtype=float
	mu, sigma = 0, 1.
	index = [[]] * Nc
	wvli = [[]] * Nc
	k = 0 ##tracing mark of wvl
	for i in range(Nc):
		fits_image_filename = 'm83_' + chnl[i] + '_cube'
		fits_ref0_filename = 'm83_' + chnl[i] + '_cube_ref'
		##zoom
		w0 = hd.rhdr(fits_image_path, fits_ref0_filename)[15]
		xmin, ymin = w0.all_world2pix(ramin, decmin, 1)
		xmax, ymax = w0.all_world2pix(ramax, decmax, 1)
		if xmin > xmax:
			t = xmax
			xmax = xmin
			xmin = t
		if ymin > ymax:
			t = ymax
			ymax = ymin
			ymin = t
		xmin = int(math.floor(xmin)) #change with chnl
		xmax = int(math.ceil(xmax))
		ymin = int(math.floor(ymin))
		ymax = int(math.ceil(ymax))
		##verify zoom validity
		if xmin < 0:
			xmin = 0
		if ymin < 0:
			ymin = 0
		if xmax < 0 or ymax < 0:
			print('zoom invalid !')
			exit(0)
		##extract data
		NAXIS1, NAXIS2, NAXIS3 = hd.rhdr3d(fits_image_path, fits_image_filename)[3:6]
		if xmax > NAXIS1:
			xmax = NAXIS1
		if ymax > NAXIS2:
			ymax = NAXIS2
		if xmin > NAXIS1 or ymin > NAXIS2:
			print('zoom invalid !')
			exit(0)
		##flux & wvl
		if i == Nc: #attention! not adapted to all (Nc-1 with akari, in convol.py)
			data0 = hd.rakari(fits_image_path, 'm83_akari_cube_raw')[0]
			data1 = hd.rakari(fits_image_path, 'm83_akari_cube_raw')[2]
		else:
			data0, data1 = hd.rdata3d(fits_image_path, fits_image_filename)
		##delete invalide slices/wvl
		if i == 1: #sl1
			del_f = 2
			del_b = 9#4
		elif i == 2: #ll2
			del_f = 4#3
			del_b = 2
		else:
			del_f = 2
			del_b = 2
		index[i] = []
		wvli[i] = []
		for jj in range(NAXIS3-del_f-del_b):
			j = jj + del_f
			dataj = data0[j, ymin:ymax, xmin:xmax]
			if i == Nc:
				if (np.isnan(dataj)).all():
					None #delete akari's wvl without data
				else:
					fits_slice_filename = 'm83_' + chnl[i] + '_' + str(j)
					hd.zoom(dataj, fits_image_path, fits_image_filename, \
						fits_slice_path, fits_slice_filename, xmin, xmax, ymin, ymax, w0)
					wvli[i].append(data1[j])
					index[i].append(j)
			else:
				fits_slice_filename = 'm83_' + chnl[i] + '_' + str(j)
				hd.zoom(dataj, fits_image_path, fits_image_filename, \
					fits_slice_path, fits_slice_filename, xmin, xmax, ymin, ymax, w0)
				wvli[i].append(data1[0][0][j][0])
				index[i].append(j)
			k += 1
	return index, wvli



def dec_pp(ramin, ramax, decmin, decmax): #dtype=float
	mu, sigma = 0, 1.
	k = 0 ##tracing mark of wvl
	for i in range(Nc):
		fits_image_filename = 'm83_' + chnl[i] + '_cube'
		fits_ref0_filename = 'm83_' + chnl[i] + '_cube_ref'
		fits_unc_filename = 'm83_' + chnl[i] + '_cube_unc'
		##zoom
		w0 = hd.rhdr(fits_image_path, fits_ref0_filename)[15]
		xmin, ymin = w0.all_world2pix(ramin, decmin, 1)
		xmax, ymax = w0.all_world2pix(ramax, decmax, 1)
		if xmin > xmax:
			t = xmax
			xmax = xmin
			xmin = t
		if ymin > ymax:
			t = ymax
			ymax = ymin
			ymin = t
		xmin = int(math.floor(xmin)) #change with chnl
		xmax = int(math.ceil(xmax))
		ymin = int(math.floor(ymin))
		ymax = int(math.ceil(ymax))
		##verify zoom validity
		if xmin < 0:
			xmin = 0
		if ymin < 0:
			ymin = 0
		##extract data
		NAXIS1, NAXIS2, NAXIS3 = hd.rhdr3d(fits_image_path, fits_image_filename)[3:6]
		if xmax > NAXIS1:
			xmax = NAXIS1
		if ymax > NAXIS2:
			ymax = NAXIS2
		NAXISx = xmax - xmin
		NAXISy = ymax - ymin
		data0 = hd.rdata(fits_image_path, fits_image_filename)
		unc = hd.rdata(fits_image_path, fits_unc_filename)
		#print(chnl[i], np.array(data).shape, np.array(unc).shape)
		##delete invalide slices/wvl
		if i == 1: #sl1
			del_f = 2
			del_b = 9#4
		elif i == 2: #ll2
			del_f = 4#3
			del_b = 2
		else:
			del_f = 2
			del_b = 2
		for jj in range(NAXIS3-del_f-del_b):
			j = jj + del_f
			dataj = data0[j, ymin:ymax, xmin:xmax]
			if i == Nc:
				if (np.isnan(dataj)).all():
					None #delete akari's wvl without data
				else:
					uncj = unc[j, ymin:ymax, xmin:xmax]
					s = np.random.normal(mu, sigma, NAXISx*NAXISy) #theta ~ N(0, 1)
					theta = s.reshape((NAXISy, NAXISx))
					dataj += theta * uncj
					fits_slice_filename = 'm83_' + chnl[i] + '_' + str(j)
					hd.zoom(dataj, fits_image_path, fits_image_filename, \
						fits_slice_path, fits_slice_filename, xmin, xmax, ymin, ymax, w0)
			else:
				uncj = unc[j, ymin:ymax, xmin:xmax]
				s = np.random.normal(mu, sigma, NAXISx*NAXISy) #theta ~ N(0, 1)
				theta = s.reshape((NAXISy, NAXISx))
				dataj += theta * uncj
				fits_slice_filename = 'm83_' + chnl[i] + '_' + str(j)
				hd.zoom(dataj, fits_image_path, fits_image_filename, \
					fits_slice_path, fits_slice_filename, xmin, xmax, ymin, ymax, w0)
			k += 1
