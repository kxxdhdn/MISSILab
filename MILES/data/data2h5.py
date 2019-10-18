#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This table was calculated by B.T. Draine, Princeton University,
Wed Dec 17 12:07:38 EST 2003
Later saved in HDF5 format by D. Hu, CEA,
Tue Aug 27 2019


Tabulated quantities:
lambda  = wavelength in vacuo (micron)
albedo  = (scattering cross section)/(extinction cross section)
<cos>   = <cos(theta)> for scattered light
C_ext/H = extinction cross section per H nucleon (cm^2/H)
K_abs   = absorption cross section per mass of dust (cm^2/gram)
<cos^2> = <cos^(theta)> for scattered light

1.870E-26 = M_dust per H nucleon (gram/H) for this dust model
1.236E+02 = M_gas/M_dust for this dust model (assuming He/H=0.096)

  lambda    albedo   <cos>  C_ext/H    K_abs   <cos^2>
 (micron)                   (cm^2/H)  (cm^2/g)          comment
----------- ------  ------ --------- --------- ------- --------

"""

from sinout import write_hdf5
import numpy as np

path = '/Users/dhu/Dropbox/MILES/data/'
fdata = path + 'ext_curve_draine03.txt'

dat = []
with open(fdata, 'r') as f:
	for line in f:
		dat.append(line)
data = np.genfromtxt(dat, dtype={'names': ('lambda', 'albedo', \
								 'cos', 'C_ext', \
								 'K_abs', 'costimes2'), \
								 'formats': ('f8', 'f8', 'f8', \
								 'f8', 'f8', 'f8')}, \
								 usecols=(0,1,2,3,4,5))

hfile = path + 'ext_curve'
write_hdf5(data['lambda'], name='Wavelength (micron)', file=hfile)
write_hdf5(data['C_ext'], name='CextovH (cm2ovH)', file=hfile, append=True)

print(">> Coucou, data2h5 done. <<")
