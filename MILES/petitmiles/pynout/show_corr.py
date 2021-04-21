#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This is the visualization of correlations (sim + fit)

"""

import os, pathlib
import numpy as np
import matplotlib.pyplot as plt

## astylo
from astylo.iolib import read_hdf5

def non_corr_err(x, y, xerr, yerr):
    '''
    f(x/y) = exp(x) / exp(y)
    df = exp(x)/exp(y) * sqrt(dx^2 + dy^2) (non-correlated)
    '''
    return np.sqrt(xerr**2 + yerr**2) * np.exp(x)/np.exp(y)
    
## Path
mroot = os.path.dirname(os.path.abspath(__file__))+'/../'
path_out = mroot+'out1/'
path_fig = path_out+'Figures/'
if not os.path.exists(path_fig):
    os.makedirs(path_fig)
filobs = path_out+'galspec'

mode = ['chi2','HB']
for m in mode:
    filfit = path_out+'fitpar_'+m
    
    fitpar = pathlib.Path(filfit+'.h5')
    if fitpar.exists():
        ## Correlations
        corrname = read_hdf5(filobs, 'Correlation name')
        Ncorr = np.size(corrname)
        
        ## Simulation
        Isim = read_hdf5(filobs, 'Correlated band intensity (Wovm2ovsr)')
        indpar = read_hdf5(filobs, 'Correlated band indpar')
        Rname = read_hdf5(filobs, 'Band ratio name')
        Nbandr, Nq, Np = Isim.shape
        Nrat = Nbandr - 1
        Rsim = np.ones((Nrat,Nq,Np))
        Rsim[0,:,:] = Isim[0,:,:]/Isim[3,:,:]
        Rsim[1,:,:] = Isim[1,:,:]/Isim[3,:,:]
        Rsim[2,:,:] = Isim[2,:,:]/Isim[3,:,:]
        Rsim[3,:,:] = Isim[4,:,:]/Isim[1,:,:]
        
        ## Chi2 (suppose all par non-correlated)
        if m=='chi2':
            par = read_hdf5(filfit, 'Best fitted parameter value')
            parerr = read_hdf5(filfit, 'Best fitted parameter error')
            lnIfit = []
            dlnIfit = []
            for i in range(Nbandr):
                lnIfit.append(par[indpar[i]-1,:,:])
                dlnIfit.append(parerr[indpar[i]-1,:,:])
            lnIfit = np.array(lnIfit)
            dlnIfit = np.array(dlnIfit)
            
            ## Band ratios
            Rfit = np.ones((Nrat,Nq,Np))
            Rerr = np.ones((Nrat,Nq,Np))
            Rfit[0,:,:] = np.exp(lnIfit[0,:,:])/np.exp(lnIfit[3,:,:])
            Rerr[0,:,:] = non_corr_err(lnIfit[0,:,:],lnIfit[3,:,:],dlnIfit[0,:,:],dlnIfit[3,:,:])
            Rfit[1,:,:] = np.exp(lnIfit[1,:,:])/np.exp(lnIfit[3,:,:])
            Rerr[1,:,:] = non_corr_err(lnIfit[1,:,:],lnIfit[3,:,:],dlnIfit[1,:,:],dlnIfit[3,:,:])
            Rfit[2,:,:] = np.exp(lnIfit[2,:,:])/np.exp(lnIfit[3,:,:])
            Rerr[2,:,:] = non_corr_err(lnIfit[2,:,:],lnIfit[3,:,:],dlnIfit[2,:,:],dlnIfit[3,:,:])
            Rfit[3,:,:] = np.exp(lnIfit[4,:,:])/np.exp(lnIfit[1,:,:])
            Rerr[3,:,:] = non_corr_err(lnIfit[4,:,:],lnIfit[1,:,:],dlnIfit[4,:,:],dlnIfit[1,:,:])
        elif m=='HB':
            Rfit = read_hdf5(filfit, 'Mean of band ratio value')
            Rerr = read_hdf5(filfit, 'Sigma of band ratio value')
        
        mlist = ['.','^','*',
                 '.','^','*',
                 '.','^','*']
        clist1 = ['lightgrey','lightgrey','lightgrey',
                  'grey','grey','grey',
                  'k','k','k']
        clist2 = ['r','r','r',
                  'g','g','g',
                  'b','b','b']
        labelist = ['c.1_SN1','c1_SN1','c10_SN1',
                    'c.1_SN25','c1_SN25','c10_SN25',
                    'c.1_SN125','c1_SN125','c10_SN125',]
        
        for i in range(Ncorr):
            # title = 'corr_'+m+'_'+str(i+1)
            title = 'corr'+str(i+1)+'_'+m
            
            filename = path_fig+title+'.png'
        
            if i==0:
                x = 0
                y = 1
            elif i==1:
                x = 0
                y = 2
            elif i==2:
                x = 0
                y = 3
            elif i==3:
                x = 1
                y = 2
            elif i==4:
                x = 1
                y = 3
            elif i==5:
                x = 2
                y = 3
                
            plt.figure(figsize=(10,6))
        
            for q in range(Nq):
                plt.scatter(Rsim[x,q,:], Rsim[y,q,:],
                             c=clist1[q], marker=mlist[q], label='sim_'+labelist[q])
                if m=='chi2':
                    plt.errorbar(Rfit[x,q,:], Rfit[y,q,:], #yerr=Rerr[y,q,:], xerr=Rerr[x,q,:],
                                c=clist2[q], fmt=mlist[q], label='fit_'+labelist[q])
                elif m=='HB':
                    plt.errorbar(Rfit[x,q,:], Rfit[y,q,:], yerr=Rerr[y,q,:], xerr=Rerr[x,q,:],
                                c=clist2[q], fmt=mlist[q], label='fit_'+labelist[q])
        
            plt.xlim(1.e-1,5.e0)
            plt.ylim(1.e-1,5.e0)
            plt.xscale('log')
            plt.yscale('log')
            plt.xlabel(Rname[x])
            plt.ylabel(Rname[y])
            plt.legend(loc='upper left')
            
            plt.savefig(filename)

# plt.show()

print('>>> Coucou show_corr [Done] <<<')
