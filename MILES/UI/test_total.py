#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This is the visualization of the test(_specModel) of auxil.f90

"""
import os
import numpy as np
import matplotlib.pyplot as plt

## astylo
from astylo.iolib import read_hdf5
from astylo.plib import plot2d_m, plot2d

## Path
mroot = os.path.dirname(os.path.abspath(__file__))+'/../'
path_out = mroot+'tests/out/'
filename = 'test_total'

## Read h5 file
wvl, FnuOBS, FnuMOD, FnuBB, FnuLINE, FnuBAND, FnuSTAR = \
    read_hdf5(path_out+filename, 
    'Wavelength (microns)', 'FnuOBS (MJyovsr)', 'FnuMOD (MJyovsr)', 
    'FnuBB (MJyovsr)', 'FnuLINE (MJyovsr)', 'FnuBAND (MJyovsr)', 
    'FnuSTAR (MJyovsr)')

p = plot2d_m(wvl, [FnuOBS, FnuMOD, FnuBB, FnuLINE, FnuBAND, FnuSTAR], 
	         xall=wvl, xlog=1, ylog=1, xlim=(2., 40.), #ylim=(0,5), 
             xlab=r'$\lambda\,(micron)$', ylab=r'$F_{\nu} \,(MJy/sr)$', 
             lablist=['Obs', 'Total', 'cont', 'lines', 'bands', 'star'], 
             legend='best', cl=['y', 'k', 'pink', 'g', 'b', 'orange'])
p.save(path_out+filename+'.png')
p.show()

print('>>> Coucou test_specModel [Done] <<<')

