#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This is the visualization of correlations (sim + fit)

"""

import os, pathlib
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

## astylo
from astylo.iolib import read_hdf5

def func(x, a):
    '''
    f(x) = a * x
    '''
    return a * x

def non_corr_df1(x, xerr):
    '''
    f(x) = exp(x)
    df = exp(x) * dx
    '''
    return xerr * np.exp(x)

def non_corr_df2(x, y, xerr, yerr):
    '''
    f(x/y) = exp(x-y)
    df = exp(x)/exp(y) * sqrt(dx^2 + dy^2) (non-correlated)
    '''
    return np.sqrt(xerr**2 + yerr**2) * np.exp(x-y)


## Path
mroot = os.path.dirname(os.path.abspath(__file__))+'/../out1/'
path_fig = mroot+'Figures/'
if not os.path.exists(path_fig):
    os.makedirs(path_fig)
filobs = mroot+'galspec'

## Correlations
corrname = read_hdf5(filobs, 'Correlation name')
Ncorr = np.size(corrname)

## Simulation
parsim = read_hdf5(filobs, 'Simulated parameter value')
indpar = read_hdf5(filobs, 'Correlated band indpar') - 1
Rname = read_hdf5(filobs, 'Band ratio name')
Nq, Np = parsim.shape[1:3]
Nbandr = len(indpar)
Nrat = len(Rname)
Rsim = np.ones((Nrat,Nq,Np))
for i in range(Nrat):
    if i==3:
        Rsim[i,:,:] = np.exp(parsim[indpar[i+1]])
    else:
        Rsim[i,:,:] = np.exp(parsim[indpar[i]])
                
mode = ['chi2','BB', 'HB']
for m in mode:
    filfit = mroot+'fitpar_'+m
    
    fitpar = pathlib.Path(filfit+'.h5')
    if fitpar.exists():
        
        Rfit = np.ones((Nrat,Nq,Np))
        Rerr = np.ones((Nrat,Nq,Np))
        
        ## Chi2 (suppose all par non-correlated)
        if m=='chi2':
            par = read_hdf5(filfit, 'Best fitted parameter value')
            parerr = read_hdf5(filfit, 'Best fitted parameter error')
            ## Band ratios
            for i in range(Nrat):
                if i==3:
                    Rfit[i,:,:] = np.exp(par[indpar[i+1],:,:])
                    Rerr[i,:,:] = non_corr_df1(par[indpar[i+1],:,:], parerr[indpar[i+1],:,:])
                else:
                    Rfit[i,:,:] = np.exp(par[indpar[i],:,:])
                    Rerr[i,:,:] = non_corr_df1(par[indpar[i],:,:], parerr[indpar[i],:,:])
            Rfit_chi2 = Rfit
            Rerr_chi2 = Rerr
        else:
            filog = mroot+'parlog_fitpar_'+m
            par = read_hdf5(filog, 'Parameter values')
            Nmcmc = read_hdf5(filog, 'Length of MCMC')[0]
            t_end = Nmcmc
            t_burnin = int(t_end/10) - 1
            ## Band ratios
            for i in range(Nrat):
                if i==3:
                    Rfit[i,:,:] = np.mean( np.exp(par[t_burnin:,indpar[i+1],:,:]), axis=0 )
                    Rerr[i,:,:] = np.std( np.exp(par[t_burnin:,indpar[i+1],:,:]), axis=0 )
                else:
                    Rfit[i,:,:] = np.mean( np.exp(par[t_burnin:,indpar[i],:,:]), axis=0 )
                    Rerr[i,:,:] = np.std( np.exp(par[t_burnin:,indpar[i],:,:]), axis=0 )
            if m=='BB':
                Rfit_BB = Rfit
                Rerr_BB = Rerr
            elif m=='HB':
                Rfit_HB = Rfit
                Rerr_HB = Rerr

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

            ## In one figure
            ##---------------
            # title = 'corr'+str(i+1)+'_all_'+m
            # filename = path_fig+title+'.png'

            # plt.figure(figsize=(10,6))
        
            # for q in range(Nq):
            #     plt.scatter(Rsim[x,q,:], Rsim[y,q,:],
            #                 c=clist1[q], marker=mlist[q], label='sim_'+labelist[q])
            #     if m=='chi2':
            #         plt.errorbar(Rfit[x,q,:], Rfit[y,q,:], yerr=Rerr[y,q,:], xerr=Rerr[x,q,:],
            #                      c=clist2[q], fmt=mlist[q], label='fit_'+labelist[q])
            #     else:
            #         plt.errorbar(Rfit[x,q,:], Rfit[y,q,:], yerr=Rerr[y,q,:], xerr=Rerr[x,q,:],
            #                      c=clist2[q], fmt=mlist[q], label='fit_'+labelist[q])
        
            # plt.xlim(1.e-1,1.e1)
            # plt.ylim(1.e-1,1.e1)
            # plt.xscale('log')
            # plt.yscale('log')
            # plt.xlabel(Rname[x])
            # plt.ylabel(Rname[y])
            # plt.legend(loc='upper left')
            
            # plt.savefig(filename)

## Compare SN 
##------------
for i in range(Ncorr):
    
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

    title1 = 'corr'+str(i+1)
    for j in range(int(Nq/3)):
        if j==0:
            title = title1+'_SN_c.1'
        elif j==1:
            title = title1+'_SN_c1'
        elif j==2:
            title = title1+'_SN_c10'
        filename = path_fig+title+'.png'
        
        fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(9,9))
        plt.subplots_adjust(left=.1, bottom=.05, \
                            right=.99, top=.95, wspace=.3, hspace=.4)
        
        for q in range(Nq):
            px, py = q%3, int(q/3)
            iq = py*3+j
            axes[px,py].scatter(Rsim[x,iq,:], Rsim[y,iq,:],
                                c=clist1[iq], marker=mlist[iq])
            ## Log linear fit
            popt, pcov = curve_fit(func, Rsim[x,iq,:], Rsim[y,iq,:])
            axes[px,py].plot(Rsim[x,iq,:], func(Rsim[x,iq,:], *popt),
                             'k-', label='a(sim)={:.2}'.format(popt[0]))
            
            if px==0:
                fitpar = pathlib.Path(mroot+'fitpar_chi2.h5')
                if fitpar.exists():
                    axes[px,py].errorbar(Rfit_chi2[x,iq,:], Rfit_chi2[y,iq,:],
                                         yerr=Rerr_chi2[y,iq,:], xerr=Rerr_chi2[x,iq,:],
                                         c=clist2[iq], fmt=mlist[iq])
                    axes[px,py].set_title('Chi2_'+labelist[iq])
                ## Log linear fit
                popt, pcov = curve_fit(func, Rfit_chi2[x,iq,:], Rfit_chi2[y,iq,:])
                axes[px,py].plot(Rfit_chi2[x,iq,:], func(Rfit_chi2[x,iq,:], *popt),
                                 'y-', label='a(fit)={:.2}'.format(popt[0]))

            elif px==1:
                fitpar = pathlib.Path(mroot+'fitpar_BB.h5')
                if fitpar.exists():
                    axes[px,py].errorbar(Rfit_BB[x,iq,:], Rfit_BB[y,iq,:],
                                         yerr=Rerr_BB[y,iq,:], xerr=Rerr_BB[x,iq,:],
                                         c=clist2[iq], fmt=mlist[iq])
                    axes[px,py].set_title('BB_'+labelist[iq])
                ## Log linear fit
                popt, pcov = curve_fit(func, Rfit_BB[x,iq,:], Rfit_BB[y,iq,:])
                axes[px,py].plot(Rfit_BB[x,iq,:], func(Rfit_BB[x,iq,:], *popt),
                                 'y-', label='a(fit)={:.2}'.format(popt[0]))
                
            elif px==2:
                fitpar = pathlib.Path(mroot+'fitpar_HB.h5')
                if fitpar.exists():
                    axes[px,py].errorbar(Rfit_HB[x,iq,:], Rfit_HB[y,iq,:],
                                         yerr=Rerr_HB[y,iq,:], xerr=Rerr_HB[x,iq,:],
                                         c=clist2[iq], fmt=mlist[iq])
                    axes[px,py].set_title('HB_'+labelist[iq])
                ## Log linear fit
                popt, pcov = curve_fit(func, Rfit_HB[x,iq,:], Rfit_HB[y,iq,:])
                axes[px,py].plot(Rfit_HB[x,iq,:], func(Rfit_HB[x,iq,:], *popt),
                                 'y-', label='a(fit)={:.2}'.format(popt[0]))
        
            axes[px,py].set_xlim(1.e-1,1.e1)
            axes[px,py].set_ylim(1.e-1,1.e1)
            axes[px,py].set_xscale('log')
            axes[px,py].set_yscale('log')
            axes[px,py].set_xlabel(Rname[x])
            axes[px,py].set_ylabel(Rname[y])
            axes[px,py].legend(loc='upper left')
        
        plt.savefig(filename)

## Compare cont
##--------------
for i in range(Ncorr):
    
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

    title1 = 'corr'+str(i+1)
    for j in range(int(Nq/3)):
        if j==0:
            title = title1+'_c_SN1'
        elif j==1:
            title = title1+'_c_SN25'
        elif j==2:
            title = title1+'_c_SN125'
        filename = path_fig+title+'.png'
        
        fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(9,9))
        plt.subplots_adjust(left=.1, bottom=.05, \
                            right=.99, top=.95, wspace=.3, hspace=.4)
        
        for q in range(Nq):
            py, px = q%3, int(q/3)
            iq = py*3+j
            axes[px,py].scatter(Rsim[x,iq,:], Rsim[y,iq,:],
                                c=clist1[iq], marker=mlist[iq])
            ## Log linear fit
            popt, pcov = curve_fit(func, Rsim[x,iq,:], Rsim[y,iq,:])
            axes[px,py].plot(Rsim[x,iq,:], func(Rsim[x,iq,:], *popt),
                             'k-', label='a(sim)={:.2}'.format(popt[0]))
            
            if px==0:
                fitpar = pathlib.Path(mroot+'fitpar_chi2.h5')
                if fitpar.exists():
                    axes[px,py].errorbar(Rfit_chi2[x,iq,:], Rfit_chi2[y,iq,:],
                                         yerr=Rerr_chi2[y,iq,:], xerr=Rerr_chi2[x,iq,:],
                                         c=clist2[iq], fmt=mlist[iq])
                    axes[px,py].set_title('Chi2_'+labelist[iq])
                ## Log linear fit
                popt, pcov = curve_fit(func, Rfit_chi2[x,iq,:], Rfit_chi2[y,iq,:])
                axes[px,py].plot(Rfit_chi2[x,iq,:], func(Rfit_chi2[x,iq,:], *popt),
                                 'y-', label='a(fit)={:.2}'.format(popt[0]))
                
            elif px==1:
                fitpar = pathlib.Path(mroot+'fitpar_BB.h5')
                if fitpar.exists():
                    axes[px,py].errorbar(Rfit_BB[x,iq,:], Rfit_BB[y,iq,:],
                                         yerr=Rerr_BB[y,iq,:], xerr=Rerr_BB[x,iq,:],
                                         c=clist2[iq], fmt=mlist[iq])
                    axes[px,py].set_title('BB_'+labelist[iq])
                ## Log linear fit
                popt, pcov = curve_fit(func, Rfit_BB[x,iq,:], Rfit_BB[y,iq,:])
                axes[px,py].plot(Rfit_BB[x,iq,:], func(Rfit_BB[x,iq,:], *popt),
                                 'y-', label='a(fit)={:.2}'.format(popt[0]))
                
            elif px==2:
                fitpar = pathlib.Path(mroot+'fitpar_HB.h5')
                if fitpar.exists():
                    axes[px,py].errorbar(Rfit_HB[x,iq,:], Rfit_HB[y,iq,:],
                                         yerr=Rerr_HB[y,iq,:], xerr=Rerr_HB[x,iq,:],
                                         c=clist2[iq], fmt=mlist[iq])
                    axes[px,py].set_title('HB_'+labelist[iq])
                ## Log linear fit
                popt, pcov = curve_fit(func, Rfit_HB[x,iq,:], Rfit_HB[y,iq,:])
                axes[px,py].plot(Rfit_HB[x,iq,:], func(Rfit_HB[x,iq,:], *popt),
                                 'y-', label='a(fit)={:.2}'.format(popt[0]))
        
            axes[px,py].set_xlim(1.e-1,1.e1)
            axes[px,py].set_ylim(1.e-1,1.e1)
            axes[px,py].set_xscale('log')
            axes[px,py].set_yscale('log')
            axes[px,py].set_xlabel(Rname[x])
            axes[px,py].set_ylabel(Rname[y])
            axes[px,py].legend(loc='upper left')
        
        plt.savefig(filename)
        
# plt.show()

print('>>> Coucou show_corr [Done] <<<')
