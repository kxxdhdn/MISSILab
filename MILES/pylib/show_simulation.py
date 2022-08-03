#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This is the visualization of simulation_MIR.h5 (simulated spectra)

"""
import os
import numpy as np
import matplotlib.pyplot as plt

## laputan
from laputan.inout import read_hdf5
from laputan.plots import pplot

## local
from librarian import croot

col_tab = ['darkred','r','pink',
           # 'orange','gold','y',
           'darkgreen','g','lime',
           'darkblue','b','cyan']

labelist = ['c1_SN5','c3_SN5','c9_SN5',
            'c1_SN20','c3_SN20','c9_SN20',
            'c1_SN80','c3_SN80','c9_SN80',]

## Path
##------
path_out = croot+'/../out/'
path_fig = path_out+'Figures/'
if not os.path.exists(path_fig):
    os.makedirs(path_fig)
filout = path_out+'simulation_MIR'

## Read h5 file
spec_unit = read_hdf5(filout, 'spectral unit')[0]
wvl = read_hdf5(filout, 'wavelength (microns)')
FnuOBS = read_hdf5(filout, 'FnuOBS ('+spec_unit+')')
dFnuOBS = read_hdf5(filout, 'dFnuOBS ('+spec_unit+')')
FnuBAND = read_hdf5(filout, 'FnuBAND ('+spec_unit+')')
Nw, Ny, Nx = FnuOBS.shape
Nband = FnuBAND.shape[0]

## Make plots
xlab = r'$\lambda\,(\mu m)$'
if spec_unit=='MKS':
    ylab = r'Surface brightness $\,(W/m2/Hz/sr)$'
elif spec_unit=='MJyovsr':
    ylab=r'Surface brightness $\,(MJy/sr)$'

## Create data tables
lab_tab = np.empty((Nx,Ny), dtype=('<U30'))
titlist = np.empty((Nx,Ny), dtype=('<U60'))
Fnu_tab = []
dFnu_tab = []
for x in range(Nx):
    for y in range(Ny):
        lab_tab[x,y] = '('+str(x+1)+','+str(y+1)+') - ['+labelist[y]+']'
        titlist[x,y] = 'Simulated Spectra '+lab_tab[x,y]+'] (rest frame)'
        Fnu_tab.append(FnuOBS[:,y,x])
        dFnu_tab.append(dFnuOBS[:,y,x])
        
Fnu_tab = np.array(Fnu_tab).reshape((Nx,Ny,Nw))
dFnu_tab = np.array(dFnu_tab).reshape((Nx,Ny,Nw))

## Plot individual spectrum with S/N
##------------------------------------
for y in range(Ny):
    for x in range(Nx):
        filename1 = path_fig+'simulation_('+str(x+1)+','+str(y+1)+').png'

        p = pplot(wvl, Fnu_tab[x,y,:]/dFnu_tab[x,y,:], yerr=0.,
                  xlog=0, ylog=0, xlim=(4.5, 22.), ylim=(0,1.e3), 
                  title=titlist[x,y], xlab=xlab, ylab=ylab+'/SN ratio', 
                  clib=['y', 'k', 'r'], label='S/N',
                  legend='best', figsize=(10,6))
        p.add_plot(wvl, Fnu_tab[x,y,:], yerr=dFnu_tab[x,y,:], label='Data')
        p.add_plot(wvl, dFnu_tab[x,y,:], yerr=0., label='Unc')
        
        p.save(filename1)

## Plot stacked spectra
##----------------------
for x in range(Nx):
    filename2 = path_fig+'simulation_stack_x='+str(x+1)+'.png'

    p = pplot(xlog=0, ylog=0, xlim=(4.5, 22.), ylim=(0.,1.e3), 
              title=titlist[x,0], xlab=xlab, ylab=ylab, clib=col_tab,
              legend='best', figsize=(10,6))
    for y in range(Ny):
        p.add_plot(wvl, Fnu_tab[x,y,:], label=lab_tab[x,y])
    
    p.save(filename2)

## Plot individual band in table
##-------------------------------
# for x in range(Nx):
#     filename3 = path_fig+'simulation_table_x='+str(x+1)+'.png'

#     fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(9,9))
#     plt.subplots_adjust(left=.1, bottom=.05, \
#                         right=.99, top=.95, wspace=.3, hspace=.4)
#     for y in range(Ny):
#         px, py = int(y/3), y%3
#         axes[px,py].plot(wvl, Fnu_tab[x,y,:],c='k')
#         for i in range(2):
#             i += 4
#             axes[px,py].plot(wvl, FnuBAND[i,:,y,x],c='b')
#         axes[px,py].set_xlim((6., 6.4))
#         axes[px,py].set_ylim((0, 2.e2))
#         axes[px,py].set_xlabel(xlab)
#         axes[px,py].set_ylabel(ylab)
#         axes[px,py].set_title(
#             fr'$Spectra ({x+1},{y+1})$')

#     plt.savefig(filename3)
    
# plt.show()

print('>>> Coucou show_simulation [Done] <<<')
