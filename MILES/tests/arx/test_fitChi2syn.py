#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This is the visualization of test_fitChi2syn

"""
import os
import numpy as np
import matplotlib.pyplot as plt

## astylo
from astylo.iolib import read_hdf5
from astylo.plib import plot2d_m

## Path
mroot = os.path.dirname(os.path.abspath(__file__))+'/'
path_out = mroot+'out/'
filename = 'test_fitChi2syn'

## Read h5 file
wvl = read_hdf5(path_out+filename, 'Wavelength (microns)')
FnuOBS = read_hdf5(path_out+filename, 'FnuOBS (MKS)')
FnuMOD = read_hdf5(path_out+filename, 'FnuMOD (MKS)')
FnuCONT = read_hdf5(path_out+filename, 'FnuCONT (MKS)')
FnuLINE = read_hdf5(path_out+filename, 'FnuLINE (MKS)')
FnuBAND = read_hdf5(path_out+filename, 'FnuBAND (MKS)')
FnuSTAR = read_hdf5(path_out+filename, 'FnuSTAR (MKS)')

p = plot2d_m(wvl, [FnuOBS, FnuMOD, FnuCONT, FnuLINE, FnuBAND, FnuSTAR], 
             xall=wvl, xlog=1, ylog=1, xlim=(2., 40.), #ylim=(0,2.), 
             xlab=r'$\lambda\,(microns)$', ylab=r'$B_{\nu} \,(W/m2/Hz/sr)$', 
             lablist=['Obs', 'Total', 'cont', 'lines', 'bands', 'star'], 
             legend='best', cl=['y', 'k', 'pink', 'g', 'b', 'orange'])
p.save(path_out+filename+'.png')

plt.show()

print('>>> Coucou test_fitChi2syn [Done] <<<')

