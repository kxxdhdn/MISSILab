#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Check if lines/features are detected by eye

"""

import os, pathlib
import numpy as np
from matplotlib.ticker import ScalarFormatter, NullFormatter

## rapyuta
from rapyuta.inout import read_hdf5
from rapyuta.plots import plotool

## local
from auxil import croot, TABLine, TABand


## Path
##------
path_out = croot+'/../out/'
path_fig = path_out+'Figures/'
if not os.path.exists(path_fig):
    os.makedirs(path_fig)
filobs = path_out+'observation_MIR'

## Read h5 file
spec_unit = read_hdf5(filobs, 'spectral unit')[0]
# spec_name = read_hdf5(filobs, 'spectrum labels')
wvl = read_hdf5(filobs, 'wavelength (microns)')
FnuOBS = read_hdf5(filobs, 'FnuOBS ('+spec_unit+')')
dFnuOBS = read_hdf5(filobs, 'dFnuOBS ('+spec_unit+')')
mask = read_hdf5(filobs, 'NaN mask')
# instr = read_hdf5(filobs, 'spectroscopic module labels')
Nw, Ny, Nx = FnuOBS.shape

## Negative to NaN
# FnuOBS[FnuOBS<0] = np.nan
# for x in range(Nx):
#     for y in range(Ny):
#         for k in range (Nw):
#             if FnuOBS[k,y,x]<0:
#                 FnuOBS[k,y,x] = -FnuOBS[k,y,x]

xlabel = r'$\lambda\,(\mu m)$'
if spec_unit=='MKS':
    ylabel = r'$F_{\nu}\ \,(W/m2/Hz/sr)$'
elif spec_unit=='MJyovsr':
    ylabel = r'$F_{\nu}\ \,(MJy/sr)$'

for x in range(Nx):
    for y in range(Ny):
        p = plotool()
        p.figure((20,8),nrows=2)
        p.set_fig(left=.1, right=.95, bottom=.1, top=.95,
                  wspace=None, hspace=.5, title=None)
        
        ## (px,py) = (icol,irow)
        px, py = 1, 1
        p.set_ax(subpos=(py,px), xlog=1, ylog=1,
                 xlim=(2.5,21), ylim=(None,None), 
                 xlabel=xlabel, ylabel=ylabel,
                 xtickfsize=20, ytickfsize=20, xfsize=20, yfsize=20,
                 title='Lines', tfsize=20)
        p.plot(wvl, FnuOBS[:,y,x],# yerr=dFnuOBS,
               c='y',
               subpos=(py,px))
        ## Detected lines
        exclines = [1,3,5,6,7,8,9,11,12,13,14,15,16,17,18,19,20,21,24,25,27,29,30,31,34,36]
        for i in range(len(TABLine)):
            if i in exclines or TABLine[i]['wave']>20:
                p.ax.axvline(x=TABLine[i]['wave'], color='grey')
            else:
                print(i, TABLine[i]['label'], TABLine[i]['wave'])
                p.ax.axvline(x=TABLine[i]['wave'], color='g')
        xtic = [2.5, 3, 3.5, 4, 5, 6, 7, 8, 10, 12, 15, 20,]# 30, 40]
        xtic_min = np.arange(2.5, 20., .1)
        p.ax.set_xticks(xtic, minor=False) # major
        p.ax.set_xticks(xtic_min, minor=True) # minor
        p.ax.xaxis.set_major_formatter(ScalarFormatter()) # major
        p.ax.xaxis.set_minor_formatter(NullFormatter()) # minor
        p.ax.legend(loc='upper left',# bbox_to_anchor=(1,1),
                    fontsize=20, framealpha=0)
        ## (px,py) = (icol,irow)
        px, py = 1, 2
        p.set_ax(subpos=(py,px), xlog=1, ylog=1,
                 xlim=(2.5,21), ylim=(None,None), 
                 xlabel=xlabel, ylabel=ylabel,
                 xtickfsize=20, ytickfsize=20, xfsize=20, yfsize=20,
                 title='Bands', tfsize=20)
        p.plot(wvl, FnuOBS[:,y,x],# yerr=dFnuOBS,
               c='y',
               subpos=(py,px))
        ## Detected bands
        print('\n')
        exclbands = [3,10,18,22,27,28,32,33]
        for i in range(len(TABand)):
            if i in exclbands or TABand[i]['wave']>20:
                p.ax.axvline(x=TABand[i]['wave'], color='grey')
            else:
                print(i, TABand[i]['label'], TABand[i]['wave'])
                p.ax.axvline(x=TABand[i]['wave'], color='b')
        xtic = [2.5, 3, 3.5, 4, 5, 6, 7, 8, 10, 12, 15, 17, 20,]# 30, 40]
        xtic_min = np.arange(2.5, 20., .1)
        p.ax.set_xticks(xtic, minor=False) # major
        p.ax.set_xticks(xtic_min, minor=True) # minor
        p.ax.xaxis.set_major_formatter(ScalarFormatter()) # major
        p.ax.xaxis.set_minor_formatter(NullFormatter()) # minor
        p.ax.legend(loc='upper left',# bbox_to_anchor=(1,1),
                    fontsize=20, framealpha=0)
        
        filename = path_fig+'check_lines_'+str(x+1)+'_'+str(y+1)+'.png'
        p.save(filename, transparent=True)
                        
print('>>> Coucou check_lines [Done] <<<')
