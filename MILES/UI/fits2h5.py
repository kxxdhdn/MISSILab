#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Convert FITS to HDF5

"""
import os
import numpy as np

## astylo
from astylo.iolib import read_fits, write_hdf5

## Path
mroot = os.path.dirname(os.path.abspath(__file__))+'/../'
path_data = mroot+'tests/dat/'
fdata = path_data + 'M83'
func = path_data + 'M83_unc'

## Read
dset = read_fits(fdata, func)

## Write
write_hdf5(fdata, 'Wavelength (microns)', dset.wave)
write_hdf5(fdata, 'LnuOBS (MJyovsr)',dset.data, append=True)
write_hdf5(fdata, 'dLnuOBS (MJyovsr)', dset.unc, append=True)

print(">>> Coucou fits2h5 [Done] <<<")
