#! /usr/bin/env python3
#-*- coding: utf-8 -*-

# Modules
import numpy as NP
import h5py as H5
import subprocess as SP

# Fake spectrum
filt_UTF8 = ["MIPS1"]
w = NP.arange(10,40,0.1)
Fnu = (w**(-2)).reshape((len(w),1,1))

# Write the HDF5 file
f = H5.File("synthetic_photometry_input.h5","w")
filt_ASCII = [n.encode("ascii","ignore") for n in filt_UTF8]
ds_filt = f.create_dataset("Filter label",dtype="S"+str(len(filt_ASCII[0])), \
data=filt_ASCII,shape=(len(filt_UTF8),))
ds_bool = f.create_dataset("(docalib,dophot)",data=[1,1])
ds_w = f.create_dataset("Wavelength (microns)",data=w)
ds_Fnu = f.create_dataset("Flux (x.Hz-1)",data=Fnu)
f.close()

# Launch the Fortran code
SP.call("synthetic_photometry")
