#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

MAIN tests

"""
# if __name__=="__main__":
# 	import sys
# 	sys.path.append('kernels')
# 	sys.path.append('utils')
# 	print(sys.path)
# or add path in __init__.py

import numpy as np
from utils.myclass import blockPrint
from utils.myfunc import celest2deg, hprint#, tprint
from utils.rwfits import *
from utils.impro import *

ref_path = '/Users/dhu/data/mosaic/SMC/'
data_path = 'test_examples/'
out_path = 'data/convolved/'
rpj_path = 'data/reprojection/'
fit_path = 'fitting/M82/Intensities/'

data_filename = 'n66_LL1_cube'
ref_filename = 'mips024'
out_ref = '_ref_'+ref_filename

bands = ['I62', 'I77', 'I86', 'I112', 'I127', 'I170']
for i in range(6):
	for bd in bands:
		hdr = read_fits(fit_path+bd+'/'+bd+'_0', False)[1]
		cen = hdr['CRVAL1'], hdr['CRVAL2']
		crop(fit_path+bd+'/'+bd+'_'+str(i), fit_path+bd+'/'+bd+'_'+str(i), cen, (15,15))
