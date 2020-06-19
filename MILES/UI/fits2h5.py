#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Convert FITS to HDF5

"""
import os
import numpy as np

## astylo
from astylo.iolib import read_fits, write_hdf5
from astylo.mlib import closest

## Path
mroot = os.path.dirname(os.path.abspath(__file__))+'/../'
path_data = mroot+'tests/dat/'
fdata = path_data + 'M83'
func = path_data + 'M83_unc'

## Read
dset = read_fits(fdata, func)

## Truncut
ind = closest(dset.wave, 21.) - 1 # at 21 microns
dset.wave = dset.wave[:ind]
dset.data = dset.data[:ind,:,:]
dset.unc = dset.unc[:ind,:,:]

## Mask
mask = np.isnan(dset.data)*1

## Write
write_hdf5(fdata, 'Wavelength (microns)', dset.wave)
write_hdf5(fdata, 'FnuOBS (MJyovsr)',dset.data, append=True)
write_hdf5(fdata, 'dFnuOBS (MJyovsr)', dset.unc, append=True)
write_hdf5(fdata, 'Mask', mask, append=True)

print(">>> Coucou fits2h5 [Done] <<<")
