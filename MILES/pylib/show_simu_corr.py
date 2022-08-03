#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This is the visualization of correlations

"""

import os, pathlib
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

## rapyuta
from rapyuta.inout import read_hdf5, h5ext
from rapyuta.plots import pplot, plotool

## local
from auxil import croot

def func(x, a, b):
    '''
    f(x) = a * x + b
    '''
    return a * x + b

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
##------
path_out = croot+'/../out/'
path_fig = path_out+'Figures/'
if not os.path.exists(path_fig):
    os.makedirs(path_fig)
# filobs = path_out+'observation_MIR' # after input_[src]_[chi2/bb/hb].py
filsim = path_out+'simulation_MIR'

## Simulation
parsim = read_hdf5(filsim, 'Simulated parameter value')
indpar = read_hdf5(filsim, 'Correlated band indpar') - 1
Ncorr = len(indpar) - 1
Rname = read_hdf5(filsim, 'Band ratio name')
Nq, Np = parsim.shape[1:3]
Nbandr = len(indpar)
Nrat = len(Rname)
Rsim = np.ones((Nrat,Nq,Np))
for i in range(Nrat):
    if i==3:
        Rsim[i,:,:] = np.exp(parsim[indpar[i+1]])
    else:
        Rsim[i,:,:] = np.exp(parsim[indpar[i]])

## Simu spec colors
clist = ['k',
         'c',
         'g',
         'm',
         'orange']
## Point form distinguishes cont
mlist = [#'.','^',
         '.','^',
         '.','^',]
## Point color (sim) distinguishes SN
clist1 = [#'lightgrey','lightgrey',
          'grey','grey',
          'k','k',]
## Point color (fit) distinguishes fitting methods
clist2 = ['r','r',
          'g','g',
          'b','b',]
labelist = ['c04_SN2',
            'c20_SN2',
            'c04_SN100',
            'c20_SN100',]

## Read fit
mode = ['chi2','bb', 'hb']
for m in mode:
    filout = path_out+'fit_'+m

    h5_out = pathlib.Path(filout+h5ext)
    if h5_out.exists():
        Rfit = np.ones((Nrat,Nq,Np))
        Rerr = np.ones((Nrat,Nq,Np))
        
        ## Chi2 (suppose all par non-correlated)
        if m=='chi2':
            par = read_hdf5(filout, 'Best fitted parameter value')
            parerr = read_hdf5(filout, 'Best fitted parameter error')
            ## Band ratios
            for i in range(Nrat):
                Rfit[i,:,:] = np.exp(par[indpar[i],:,:])
                Rerr[i,:,:] = non_corr_df1(par[indpar[i],:,:], parerr[indpar[i],:,:])
            Rfit_chi2 = Rfit
            Rerr_chi2 = Rerr
        else:
            filmcmc = path_out+'parlog_fit_'+m
            par = read_hdf5(filmcmc, 'Parameter values')
            # Nmcmc = read_hdf5(filmcmc, 'Length of MCMC')[0]
            Nmcmc = read_hdf5(filmcmc, 'Last index')[0]
            t_end = Nmcmc
            t_burnin = int(t_end/10) - 1
            ## Band ratios
            for i in range(Nrat):
                Rfit[i,:,:] = np.mean( np.exp(par[t_burnin:t_end,indpar[i],:,:]), axis=0 )
                Rerr[i,:,:] = np.std( np.exp(par[t_burnin:t_end,indpar[i],:,:]), axis=0 )
            if m=='bb':
                Rfit_bb = Rfit
                Rerr_bb = Rerr
            elif m=='hb':
                Rfit_hb = Rfit
                Rerr_hb = Rerr

        ## In one figure
        ##---------------
        for i in range(Ncorr):
            
            if i==0:
                x = 0 # 3.3/11.2
                y = 1 # 6.2/11.2
            elif i==1:
                x = 0
                y = 4 # 12.7/11.2
            elif i==2:
                x = 2 # 7.7/11.2
                y = 1
            elif i==3:
                x = 2
                y = 3 # 8.6/11.2
            elif i==4:
                x = 2
                y = 4

            title = 'corr'+str(i+1)+'_all_'+m
            filename = path_fig+title+'.png'

            p = pplot(xlabel=Rname[x], ylabel=Rname[y],
                      xlog=1, ylog=1,
                      # xlim=(,), ylim=(,),
                      title=None,
                      figsize=(12,8), right=.7, left=.1, bottom=.1, top=.95,
                      legend='upper left', anchor=(1,1),
                      titlesize=20, labelsize=20, ticksize=20, legendsize=20)
            for q in range(Nq):
                p.add_plot(Rsim[x,q,:], Rsim[y,q,:],
                           fmt='s', c=clist1[q], marker=mlist[q], label='sim_'+labelist[q],)
                p.add_plot(Rfit[x,q,:], Rfit[y,q,:], yerr=Rerr[y,q,:], xerr=Rerr[x,q,:],
                           fmt='s', c=clist2[q], marker=mlist[q], label='fit_'+labelist[q])
                
            p.save(filename)#, transparent=True)

N1 = 2 # S/N
N2 = 2 # cont
xgrid = np.logspace(-4,3,10000)

## Compare SN 
##------------
for i in range(Ncorr):
    
    if i==0:
        x = 0 # 3.3/11.2
        y = 1 # 6.2/11.2
    elif i==1:
        x = 0
        y = 4 # 12.7/11.2
    elif i==2:
        x = 2 # 7.7/11.2
        y = 1
    elif i==3:
        x = 2
        y = 3 # 8.6/11.2
    elif i==4:
        x = 2
        y = 4

    title1 = 'corr'+str(i+1)
    for pz in range(N1):
        if pz==0:
            title = title1+'_SN_c04'
        elif pz==1:
            title = title1+'_SN_c20'
        filename = path_fig+title+'.png'

        p1 = plotool()
        p1.figure((12,12), nrows=3, ncols=N1)
        p1.set_fig(left=.1, bottom=.1, right=.8, top=.95,
                   wspace=1.5, hspace=.3)
        
        for py in range(int(Nq/N1)):
            iq = py*N1+pz
            for px in range(3):
                p1.set_ax((px+1,py+1), xlog=1, ylog=1,
                          xlabel=Rname[x], ylabel=Rname[y],
                          xtickfsize=15, ytickfsize=15,
                          xfsize=15, yfsize=15) # size
                p1.plot(Rsim[x,iq,:], Rsim[y,iq,:],
                        fmt='s', c=clist1[iq], marker=mlist[iq])

                ## Log linear fit
                popt, pcov = curve_fit(func, Rsim[x,iq,:], Rsim[y,iq,:])
                p1.plot(xgrid, func(xgrid, *popt),
                        c='k', label='a(sim)={:.2}'.format(popt[0]))

                ix = px*N1+py
                ## Chi2
                if px==0:
                    h5_out = pathlib.Path(path_out+'fit_chi2'+h5ext)
                    if h5_out.exists():
                        p1.plot(Rfit_chi2[x,iq,:], Rfit_chi2[y,iq,:],
                                yerr=Rerr_chi2[y,iq,:], xerr=Rerr_chi2[x,iq,:],
                                fmt='s', c=clist2[ix], marker=mlist[iq])
                        p1.ax.set_title('Chi2_'+labelist[iq])
                        ## Log linear fit
                        popt, pcov = curve_fit(func, Rfit_chi2[x,iq,:], Rfit_chi2[y,iq,:])
                        p1.plot(xgrid, func(xgrid, *popt),
                                c='y', label='a(fit)={:.2}'.format(popt[0]))
                ## BB
                elif px==1:
                    h5_out = pathlib.Path(path_out+'fit_bb'+h5ext)
                    if h5_out.exists():
                        p1.plot(Rfit_bb[x,iq,:], Rfit_bb[y,iq,:],
                                yerr=Rerr_bb[y,iq,:], xerr=Rerr_bb[x,iq,:],
                                fmt='s', c=clist2[ix], marker=mlist[iq])
                        p1.ax.set_title('BB_'+labelist[iq])
                        ## Log linear fit
                        popt, pcov = curve_fit(func, Rfit_bb[x,iq,:], Rfit_bb[y,iq,:])
                        p1.plot(xgrid, func(xgrid, *popt),
                                c='y', label='a(fit)={:.2}'.format(popt[0]))
                ## HB
                elif px==2:
                    h5_out = pathlib.Path(path_out+'fit_hb'+h5ext)
                    if h5_out.exists():
                        p1.plot(Rfit_hb[x,iq,:], Rfit_hb[y,iq,:],
                                yerr=Rerr_hb[y,iq,:], xerr=Rerr_hb[x,iq,:],
                                fmt='s', c=clist2[ix], marker=mlist[iq])
                        p1.ax.set_title('HB_'+labelist[iq])
                        ## Log linear fit
                        popt, pcov = curve_fit(func, Rfit_hb[x,iq,:], Rfit_hb[y,iq,:])
                        p1.plot(xgrid, func(xgrid, *popt),
                                c='y', label='a(fit)={:.2}'.format(popt[0]))
        
                p1.set_legend(loc='upper left', bbox_to_anchor=(1,1),
                              fontsize=15, framealpha=0) # size
        p1.save(filename)
exit()
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

    title2 = 'corr'+str(i+1)
    for pz in range(int(Nq/3)):
        if pz==0:
            title = title2+'_c_SN5'
        elif pz==1:
            title = title2+'_c_SN20'
        elif pz==2:
            title = title2+'_c_SN80'
        filename = path_fig+title+'.png'
        
        fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(9,9))
        plt.subplots_adjust(left=.1, bottom=.05, \
                            right=.99, top=.95, wspace=.3, hspace=.4)
        
        for py in range(3):
            iq = pz*3+py
            for px in range(3):
                axes[px,py].scatter(Rsim[x,iq,:], Rsim[y,iq,:],
                                    c=clist1[iq], marker=mlist[iq])
                ## Log linear fit
                popt, pcov = curve_fit(func, Rsim[x,iq,:], Rsim[y,iq,:])
                axes[px,py].plot(Rsim[x,iq,:], func(Rsim[x,iq,:], *popt),
                                 'k-', label='a(sim)={:.2}'.format(popt[0]))
                
                ix = px*3+py
                ## Chi2
                if px==0:
                    h5_out = pathlib.Path(path_out+'fit_chi2'+h5ext)
                    if h5_out.exists():
                        axes[px,py].errorbar(Rfit_chi2[x,iq,:], Rfit_chi2[y,iq,:],
                                             yerr=Rerr_chi2[y,iq,:], xerr=Rerr_chi2[x,iq,:],
                                             c=clist2[ix], fmt=mlist[iq])
                        axes[px,py].set_title('Chi2_'+labelist[iq])
                    ## Log linear fit
                    popt, pcov = curve_fit(func, Rfit_chi2[x,iq,:], Rfit_chi2[y,iq,:])
                    axes[px,py].plot(Rfit_chi2[x,iq,:], func(Rfit_chi2[x,iq,:], *popt),
                                     'y-', label='a(fit)={:.2}'.format(popt[0]))
                ## BB
                elif px==1:
                    h5_out = pathlib.Path(path_out+'fit_bb'+h5ext)
                    if h5_out.exists():
                        axes[px,py].errorbar(Rfit_bb[x,iq,:], Rfit_bb[y,iq,:],
                                             yerr=Rerr_bb[y,iq,:], xerr=Rerr_bb[x,iq,:],
                                             c=clist2[ix], fmt=mlist[iq])
                        axes[px,py].set_title('BB_'+labelist[iq])
                    ## Log linear fit
                    popt, pcov = curve_fit(func, Rfit_bb[x,iq,:], Rfit_bb[y,iq,:])
                    axes[px,py].plot(Rfit_bb[x,iq,:], func(Rfit_bb[x,iq,:], *popt),
                                     'y-', label='a(fit)={:.2}'.format(popt[0]))
                ## HB
                elif px==2:
                    h5_out = pathlib.Path(path_out+'fit_hb'+h5ext)
                    if h5_out.exists():
                        axes[px,py].errorbar(Rfit_hb[x,iq,:], Rfit_hb[y,iq,:],
                                             yerr=Rerr_hb[y,iq,:], xerr=Rerr_hb[x,iq,:],
                                             c=clist2[ix], fmt=mlist[iq])
                        axes[px,py].set_title('HB_'+labelist[iq])
                    ## Log linear fit
                    popt, pcov = curve_fit(func, Rfit_hb[x,iq,:], Rfit_hb[y,iq,:])
                    axes[px,py].plot(Rfit_hb[x,iq,:], func(Rfit_hb[x,iq,:], *popt),
                                     'y-', label='a(fit)={:.2}'.format(popt[0]))
            
                axes[px,py].set_xlim(1.e-2,1.e1)
                axes[px,py].set_ylim(1.e-1,1.e1)
                axes[px,py].set_xscale('log')
                axes[px,py].set_yscale('log')
                axes[px,py].set_xlabel(Rname[x])
                axes[px,py].set_ylabel(Rname[y])
                axes[px,py].legend(loc='upper left')
        
        plt.savefig(filename)
        
# plt.show()

print('>>> Coucou show_sim_corr [Done] <<<')
