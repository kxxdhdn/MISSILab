#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This is the visualization of the test_fitM83.h5

"""

import os
import numpy as np
import matplotlib.pyplot as plt

## astylo
from astylo.iolib import read_hdf5
from astylo.plib import plot2d_m, plot2d

## Path
mroot = os.path.dirname(os.path.abspath(__file__))+'/../'
path_dat = mroot+'tests/dat/'
path_out = mroot+'tests/out/'
path_fig = path_out+'Figures/'
if not os.path.exists(path_fig):
    os.makedirs(path_fig)

filename = 'M83'

## Read h5 files
wvl, FnuOBS, dFnuOBS, FnuMOD, FnuBB, FnuLINE, FnuBAND, FnuSTAR = \
    read_hdf5(path_out+'test_fit'+filename, 
        'Wavelength (microns)', 'FnuOBS (MJyovsr)', 'FnuUNC (MJyovsr)', 
        'FnuMOD (MJyovsr)', 'FnuBB (MJyovsr)', 'FnuLINE (MJyovsr)', 
        'FnuBAND (MJyovsr)', 'FnuSTAR (MJyovsr)')
Fbb = []
for i in range(4):
    Fbb.append(read_hdf5(path_out+'test_fit'+filename, 
        'FnuBB'+str(i+1)+' (MJyovsr)')[0])
Fline = []
for i in range(11):
    Fline.append(read_hdf5(path_out+'test_fit'+filename, 
        'FnuLINE'+str(i+1)+' (MJyovsr)')[0])
Fband = []
for i in range(31):
    Fband.append(read_hdf5(path_out+'test_fit'+filename,
        'FnuBAND'+str(i+1)+' (MJyovsr)')[0])

Nw, Ny, Nx = FnuOBS.shape
# x, y = 0, 2
for x in range(Nx):
    for y in range(Ny):
        if not any(np.isnan(FnuOBS[:,y,x])):
            errlist = [dFnuOBS[:,y,x]]
            for i in range(48):
                errlist.append(np.zeros(FnuOBS.shape[0]))

            ## Curves & colors
            curves = [FnuOBS[:,y,x], FnuMOD[:,y,x]]
            cl = ['k', 'grey', 'orange', 'pink', 'orange', 'orange']
            for i in range(4):
                curves.append(Fbb[i][:,y,x])

            for i in range(11):
                curves.append(Fline[i][:,y,x])
                cl.append('g')
            for i in range(31):
                curves.append(Fband[i][:,y,x])
                cl.append('b')
            curves.append(FnuSTAR[:,y,x])
            cl.append('y')
            p = plot2d_m(wvl, curves, yerr = errlist, 
            	         xall=wvl, xlog=0, ylog=0, xlim=(5.,23.), #ylim=(0,5), 
                         xlab=r'$\lambda\,(micron)$', ylab=r'$F_{\nu} \,(MJy/sr)$', 
                         # lablist=['Obs', 'Total', 'cont', 'lines', 'bands', 'star'], 
                         legend='best', cl=cl)
            p.save(path_fig+filename+'_('+str(x+1)+','+str(y+1)+')'+'.png')
# p.show()

print('>>> Coucou test_fitChi2 [Done] <<<')
