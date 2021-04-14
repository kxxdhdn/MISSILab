#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Tracking parameter distribution during Gibbs sampling

"""
import os
import numpy as np
import matplotlib.pyplot as plt

## astylo
from astylo.iolib import read_hdf5

## local
from utilities import TABand

## Path
mroot = os.path.dirname(os.path.abspath(__file__))+'/../out1/'
path_fig = mroot+'Figures/'
if not os.path.exists(path_fig):
    os.makedirs(path_fig)
filobs = mroot+'galspec'
filfit = mroot+'fitpar_BB'
filog = mroot+'parlog_fitpar_BB'

Nmcmc = read_hdf5(filog, 'Length of MCMC')[0]
t_end = Nmcmc
t_burnin = int(t_end/10) - 1
parname = read_hdf5(filog, 'Parameter label')
parmcmc = read_hdf5(filog, 'Parameter values')
Npar, Ny, Nx = parmcmc.shape[1:4]
parmu = read_hdf5(filfit, 'Mean of parameter value')
# mu = np.mean(parmcmc[t_burnin:t_end,:,:,:], axis=0)
# print(mu-parmu)
parsig = read_hdf5(filfit, 'Sigma of parameter value')
# sig = np.std(parmcmc[t_burnin:t_end,:,:,:], axis=0)
# print(sig-parsig)
# exit()

filmod = mroot+'input_fitMIR_model'
fixed = read_hdf5(filmod, 'parinfo fixed')
parini = read_hdf5(filmod, 'parinfo value')
labB = read_hdf5(filmod, 'label band')
parsim = read_hdf5(filobs, 'Simulated parameter value') # True par
chi2ini = read_hdf5(filobs, 'Chi2 fitted parameter value') # Chi2 result

            
## Parameter track (pix)
##-----------------------

## Find band param
bname = 'Main 6.2 (1)'
ipar = 0
bex = 0 # Cband: +0 WSband: +1; WLband: +2
for i, name in enumerate(parname):
    if name[:5]=='Cband':
        indB = int(name[5:])
        if labB[indB-1]==bname:
            ipar = i + bex
            
counter = np.arange(Nmcmc)
for x in range(Nx):
    filename = path_fig+'partrack_x='+str(x+1)+'.png'
    
    fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(9,9))
    plt.subplots_adjust(left=.1, bottom=.05, \
                        right=.99, top=.95, wspace=.3, hspace=.4)
    for y in range(Ny):
        px, py = int(y/3), y%3
        axes[px,py].scatter(parmcmc[:,ipar,y,x], counter, s=.5)
        for band in TABand:
            if band['label']==bname:
                if bex==0:
                    axes[px,py].vlines(band['wave'],0,Nmcmc,colors='r')
                elif bex==1:
                    axes[px,py].vlines(band['sigmaS'],0,Nmcmc,colors='r')
                elif bex==2:
                    axes[px,py].vlines(band['sigmaL'],0,Nmcmc,colors='r')
        axes[px,py].set_xlabel('Parameter value')
        axes[px,py].set_ylabel('MCMC counter')
        axes[px,py].set_title(
            fr'${parname[ipar]}({x+1},{y+1})\ [{parmu[ipar,y,x]:.2f},{parsig[ipar,y,x]:.2f}]$')
    
    plt.savefig(filename)
    
## Parameter distribution (pix)
##------------------------------
# Nbin = int(Nmcmc/10)
# for x in range(Nx):
#     filename = path_fig+'pardist_x='+str(x+1)+'.png'

#     fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(9,9))
#     plt.subplots_adjust(left=.1, bottom=.05, \
#                         right=.99, top=.95, wspace=.3, hspace=.4)
#     for y in range(Ny):
#         px, py = int(y/3), y%3
#         n, bins, patches = axes[px,py].hist(parmcmc[:,ipar,y,x], Nbin, facecolor='g', alpha=.5)
#         for band in TABand:
#             if band['label']=='Main 6.2 (1)':
#                 if bex==0:
#                     axes[px,py].vlines(band['wave'],0,Nmcmc,colors='r')
#                 elif bex==1:
#                     axes[px,py].vlines(band['sigmaS'],0,Nmcmc,colors='r')
#                 elif bex==2:
#                     axes[px,py].vlines(band['sigmaL'],0,Nmcmc,colors='r')
#         axes[px,py].set_xlabel('Parameter value')
#         axes[px,py].set_ylabel('Counts')
#         axes[px,py].set_title(
#             fr'${parname[ipar]}({x+1},{y+1})\ [{parmu[ipar,y,x]:.2f},{parsig[ipar,y,x]:.2f}]$')
    
#     plt.savefig(filename)

## Parameter track (par)
##-----------------------
x,y = 2,5
for j in range(int(Npar/9)+1):
    filename = path_fig+'partrack_ipar='+str(j*9+1)+'-'+str(j*9+9)+'.png'
    
    fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(9,9))
    plt.subplots_adjust(left=.1, bottom=.05, \
                        right=.99, top=.95, wspace=.3, hspace=.4)

    for jpar in range(9):
        jpar += j*9
        px, py = int(jpar%9/3), jpar%9%3
        if jpar<Npar:
            axes[px,py].scatter(parmcmc[:,jpar,y,x], counter, s=.5)

            ## Fitted (and after burning) param in magenta
            axes[px,py].vlines(parmu[jpar,y,x],0,Nmcmc,colors='m')
            axes[px,py].vlines(parmu[jpar,y,x]-parsig[jpar,y,x],0,Nmcmc,colors='m',linestyles='dashed')
            axes[px,py].vlines(parmu[jpar,y,x]+parsig[jpar,y,x],0,Nmcmc,colors='m',linestyles='dashed')

            ## Initial param in red
            if fixed[jpar]=='F':
                # axes[px,py].vlines(parini[jpar],0,Nmcmc,colors='r')
                axes[px,py].vlines(chi2ini[jpar,y,x],0,Nmcmc,colors='r')
                axes[px,py].vlines(parsim[jpar,y,x],0,Nmcmc,colors='g')
            
            axes[px,py].set_xlabel('Parameter value')
            axes[px,py].set_ylabel('MCMC counter')
            axes[px,py].set_title(
                fr'${parname[jpar]}({x+1},{y+1})\ [{parmu[jpar,y,x]:.2f},{parsig[jpar,y,x]:.2f}]$')
    
    plt.savefig(filename)

# plt.show()

print('>>> Coucou partrack [Done] <<<')
