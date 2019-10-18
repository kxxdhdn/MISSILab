#!/usr/bin/env python
# -*- coding: utf-8 -*-

from sinout import read_fits, write_hdf5
import numpy as np

path = '/Users/dhu/Dropbox/MILES/input/'
fdata = path + 'm83'
func = path + 'm83_unc'
Fnu, wave = read_fits(fdata)[0:2]
dFnu = read_fits(func)[0]

hfile = path + 'test_sample'
write_hdf5(wave, name='Wavelength (micron)', file=hfile)
write_hdf5(Fnu, name='LnuOBS (MJyovsr)', file=hfile, append=True)
write_hdf5(dFnu, name='dLnuOBS (MJyovsr)', file=hfile, append=True)

print(">> Coucou, fits2h5 done. <<")
