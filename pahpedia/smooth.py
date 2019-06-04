#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

PSF homogeneisation

"""

import numpy as np
import subprocess as SP
from kernels.gen_kern import kern
from utils.impro import *
from utils.rwcsv import write_csv

psf = [2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7.]

def crop3D(filIN, filOUT, uncIN, cen, size, wmod=0):

	## slice cube
	if uncIN!=None:
		wvl, hdr = cubislice(filIN, filOUT, uncIN, None, wmod)
	else:
		wvl, hdr = cubislice(filIN, filOUT, None, None, wmod)

	SList = []
	for k in range(np.size(wvl)):
		SList.append(filOUT+'_'+'0'*(4-len(str(k)))+str(k))
	## crop slices
	cube=[]
	for SLout in SList:
		cube.append(crop(SLout, SLout, cen, size))
	cube = np.array(cube)

	return cube, SList, wvl, hdr

def marker(filIN):
	
	# ker_name = 'Kernel_HiRes_IRAC_5.8_to_Gauss_06.0'
	ker_name = 'Kernel_HiRes_IRAC_8.0_to_Gauss_06.0'
	## write filename lists
	kern_path = 'kernels/'
	kern_filename = 'kern'
	write_csv(kern_path+kern_filename, \
		["Images", "Kernels"], [[filIN, ker_name], [filIN+'_unc', ker_name]])

def choker(SList, wvl):

	kerl = []
	sigma_lam = kern(wvl)[0]
	##choose kernels
	for k in range(np.size(wvl)):
		for j in range(np.size(psf)):
			dp = sigma_lam[k] - psf[j]
			dp1 = sigma_lam[k] - psf[j+1]
			if dp1 > 0:
				flag = 0
			elif dp > 0 and dp1 < 0:
				if dp + dp1 <= 0:
					flag = 1
					break
				else:
					j += 1
					flag = 1
					break
			else:
				print('Error: need smaller PSF! ')
				flag = 1
				break
		#print(psf[j])
		if (flag==0):
			print('Error: need bigger PSF! (In this case IndexError will appear first at psf[j+1])')
		# ker_name = 'Kernel_HiRes_Gauss_06.0_to_MIPS_24'
		ker_name = 'Kernel_HiRes_Gauss_0' + str(psf[j]) + '_to_MIPS_24'
		# ker_name = 'Kernel_HiRes_Gauss_0' + str(psf[j]) + '_to_Gauss_06.0'
		im_name = SList[k]
		kern_k = [im_name, ker_name]
		kerl.append(kern_k) # kerl-----im, kern
	
	## write filename lists
	kern_path = 'kernels/'
	kern_filename = 'kern'
	write_csv(kern_path+kern_filename, ["Images", "Kernels"], kerl)

def do_conv():

	SP.call('cd IDL/\nidl conv.pro', shell=True)

"""
------------------------------ MAIN (test) ------------------------------
"""
if __name__ == "__main__":

	## in situ test
	from astropy import units as u
	import matplotlib.pyplot as plt
	from utils.myfunc import celest2deg, fclean
	from utils.mapping import *

	ref_path = '/Users/dhu/data/mosaic/SMC/'
	data_path = 'data/n66/'
	out_path = 'data/n66/slices/'
	rpj_path = 'data/n66/reprojection/'

	data_filename = 'n66_LL1'
	ref_filename = 'mips024'
	out_ref = '_ref_'+ref_filename

	ra, dec = celest2deg(0., 59., 3.5623, -72., 10., 33.972)
	dx, dy = 34, 40
	
	## 3D cube cropping
	cube, SList, wvl, hdr = crop3D(data_path+data_filename, out_path+data_filename, None, \
		(ra, dec), (dx, dy), 1)
	print("Cropped cube size: ", cube.shape)

	choker(SList, wvl)

	do_conv()
	
	## mips data
	crop(ref_path+ref_filename, data_path+out_ref, \
		(ra, dec), (dx, dy))
	print("Cropped ref size: ", ref.shape)

	ref = rpj(data_path+out_ref, rpj_path+out_ref, SList[10]+'_conv', (dy, dx))[0]
	
	data, ft, w = rpj(SList[10]+'_conv', rpj_path+data_filename, SList[10]+'_conv', (dy, dx))
	imview([data, ref, ft], w, (1,3), figsize=(12,5))

	plt.show()
