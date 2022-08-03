#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This is the visualization of simulation_MIR.h5 (simulated spectra)

"""
import os
import numpy as np
import matplotlib.pyplot as plt

## rapyuta
from rapyuta.inout import read_hdf5
from rapyuta.plots import pplot

## local
from auxil import croot

clist = ['k',
         'c',
         'g',
         'm',
         'orange']

labelist = ['c04_SN2',
            'c20_SN2',
            'c04_SN100',
            'c20_SN100',]

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
lab_tab = np.empty((Ny,Nx), dtype=('<U30'))
titlist = np.empty((Ny,Nx), dtype=('<U60'))
Fnu_tab = []
dFnu_tab = []
for x in range(Nx):
    for y in range(Ny):
        lab_tab[y,x] = labelist[y]+'_'+str(x+1)
        titlist[y,x] = 'Synthetic Spectra '+lab_tab[y,x]+' (rest frame)'
        Fnu_tab.append(FnuOBS[:,y,x])
        dFnu_tab.append(dFnuOBS[:,y,x])
        
Fnu_tab = np.array(Fnu_tab).reshape((Nx,Ny,Nw))
dFnu_tab = np.array(dFnu_tab).reshape((Nx,Ny,Nw))

## Plot individual spectrum with S/N
##------------------------------------
for y in range(Ny):
    for x in range(Nx):
        filename1 = path_fig+'simu_'+str(x+1)+'_'+str(y+1)+'.png'

        p = pplot(xlog=1, ylog=1, nonposx='clip', nonposy='clip',
                  xlim=(2.4,21), ylim=(1e-1,1e4),
                  xlabel=r'${\rm Wavelengths}\ \lambda\ (\mu m)$',
                  ylabel=r'${\rm Surface\ brightness}\ F_{\nu}\ (MJy/sr)$',
                  title=titlist[y,x],
                  figsize=(10,8), right=.95, left=.1, bottom=.15,
                  # legend='upper left', anchor=(1,1),
                  titlesize=20, labelsize=20, ticksize=20, legendsize=20)
        p.add_plot(wvl, Fnu_tab[x,y,:], yerr=dFnu_tab[x,y,:],
                   ec='grey', elw=1, c='k', label='Data')
        p.add_plot(wvl, dFnu_tab[x,y,:], c='r', label='Errors')

        # p.ax.text(.9,.05,'(a)',size=20,c='grey',transform=p.ax.transAxes)
        p.ax.legend(loc='upper left',# bbox_to_anchor=(1,1),
                    fontsize=20, framealpha=0,)
        p.save(filename1, transparent=True)

## Plot stacked spectra
##----------------------
for x in range(Nx):
    filename2 = path_fig+'simu_stack_x='+str(x+1)+'.png'

    p = pplot(xlog=1, ylog=1, nonposx='clip', nonposy='clip',
              xlim=(2.4,21), ylim=(1e0,2e4),
              xlabel=r'${\rm Wavelengths}\ \lambda\ (\mu m)$',
              ylabel=r'${\rm Surface\ brightness}\ F_{\nu}\ (MJy/sr)$',
              title=None, clib=clist,
              figsize=(10,8), right=.95, left=.15, bottom=.15, top=.95,
              # legend='upper left', anchor=(1,1),
              titlesize=20, labelsize=20, ticksize=20, legendsize=20)
    for y in range(Ny):
        p.add_plot(wvl, Fnu_tab[x,y,:], label=lab_tab[y,x])

    # p.ax.text(.9,.05,'(a)',size=20,c='grey',transform=p.ax.transAxes)
    p.ax.legend(loc='upper left',# bbox_to_anchor=(1,1),
                fontsize=20, framealpha=0,)
    p.save(filename2, transparent=True)

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
