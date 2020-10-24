#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Convert ASCII to HDF5

"""
import os
import numpy as np

## astylo
from astylo.iolib import write_hdf5, ascext
ascext = '.txt'

## Path
mroot = os.path.dirname(os.path.abspath(__file__))+'/../'
datdir = mroot+'data/'
filename = datdir + 'extcurve'

## Read
dset = np.genfromtxt(filename+ascext, skip_header=80,
                     dtype={'names': ('lambda','albedo','cos','C_extovH','K_abs','cos2'),
                            'formats':  ('f8','f8','f8','f8','f8','f8')},
                     usecols=(0,1,2,3,4,5))

## Write
write_hdf5(filename, 'lambda (micron)', dset['lambda'])
write_hdf5(filename, 'albedo', dset['albedo'], append=True)
write_hdf5(filename, '<cos>', dset['cos'], append=True)
write_hdf5(filename, 'C_extovH (cm^2ovH)', dset['C_extovH'], append=True)
write_hdf5(filename, 'K_abs (cm^2ovg)', dset['K_abs'], append=True)
write_hdf5(filename, '<cos^2>', dset['cos2'], append=True)

print(">>> Coucou asc2h5 [Done] <<<")
