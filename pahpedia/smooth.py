#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

PSF homogeneisation

"""

#import sys
#sys.path.append('..')
#print(sys.path)
import numpy as np
import subprocess as SP
from kernels.gen_kern import kern
from impro import crop, cubislice
from utils.rwcsv import write_csv

psf = [1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8., 8.5, 9., 9.5, 10.]

def crop3D(filIN, filOUT, centre, size, nopr=True):

	## slice cube
	wvl = cubislice(filIN, filOUT, '_')
	SList = []
	for k in range(np.size(wvl)):
		SList.append(filOUT+'_'+'0'*(4-len(str(k)))+str(k)+'_')
	## crop slices
	cube=[]
	for slOUT in SList:
		cube.append(crop(slOUT, slOUT, centre, size, nopr))
	cube = np.array(cube)

	return cube, wvl, SList

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
		ker_name = 'Kernel_HiRes_Gauss_0' + str(psf[j]) + '_to_Gauss_06.0'
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
	from utils.myfunc import deg2world, fclean

	data_path = 'test_examples/'
	out_path = 'data/convolved/'
	data_filename = 'n66_LL1_cube'

	ra, dec = deg2world(0., 59., 3.5623, -72., 10., 33.972)
	d_ra, d_dec = 1./30., 1./30.
	
	## 3D cube cropping
	cube, wvl, SList = crop3D(data_path+data_filename, out_path+data_filename, [ra, dec], [d_ra, d_dec])
	print("Cropped cube size: ", cube.shape)

	choker(SList, wvl)

	do_conv()

	fclean('data/convolved/*_.fits')
#	fclean('data/convolved/*conv.fits')
