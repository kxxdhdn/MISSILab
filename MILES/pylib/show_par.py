#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Tracking parameter distribution during Gibbs sampling

"""

import os, pathlib
import numpy as np
import matplotlib.pyplot as plt

## laputan
from laputan.inout import read_hdf5

## local
from auxil import croot, TABand

## Path
##------
path_out = croot+'/../out/'
path_fig = path_out+'Figures/'
if not os.path.exists(path_fig):
    os.makedirs(path_fig)
filobs = path_out+'observation_MIR'
filsim = path_out+'simulation_MIR'
h5_sim = pathlib.Path(filsim+'.h5')

## Read h5 file
mode = ['bb','hb']
for m in mode:
    filout = path_out+'fit_'+m
    filmcmc = path_out+'parlog_fit_'+m

    h5_out = pathlib.Path(filout+'.h5')
    if h5_out.exists():
        # Nmcmc = read_hdf5(filmcmc, 'Length of MCMC')[0]
        # counter = np.arange(Nmcmc)
        Nmcmc = read_hdf5(filmcmc, 'Last index')[0]
        counter = np.arange(Nmcmc)
        t_end = Nmcmc
        t_burnin = int(t_end/10) - 1
        parname = read_hdf5(filmcmc, 'Parameter label')
        parmcmc = read_hdf5(filmcmc, 'Parameter values')
        Npar, Ny, Nx = parmcmc.shape[1:4]
        meanpar = read_hdf5(filout, 'Mean of parameter value')
        # mu = np.mean(parmcmc[t_burnin:t_end,:,:,:], axis=0)
        # print(mu-meanpar)
        stdevpar = read_hdf5(filout, 'Sigma of parameter value')
        # sig = np.std(parmcmc[t_burnin:t_end,:,:,:], axis=0)
        # print(sig-stdevpar)
        # exit()
        qpar = read_hdf5(filout, 'Quantiles of parameter value')
        
        filmod = path_out+'input_model'
        fixed = read_hdf5(filmod, 'parinfo fixed')
        parini = read_hdf5(filmod, 'parinfo value')
        labB = read_hdf5(filmod, 'label band')
        if h5_sim.exists():
            simpar = read_hdf5(filsim, 'Simulated parameter value') # True par
        chi2ini = read_hdf5(filobs, 'Chi2 fitted parameter value') # Chi2 result
        
        bname = 'Main 11.2'
        ipar = 0
        bex = 0 # Cband: +1 WSband: +2; WLband: +3
        for i, name in enumerate(parname):
            if name[:5]=='Cband':
                indB = int(name[5:])
                if labB[indB-1]==bname:
                    ipar = i + bex -1
        
        # ipar = 6 # lnFcont4
        
        ## Parameter track (pix)
        ##-----------------------
        for x in range(Nx):
            filename = path_fig+'partrack_'+m+'_x='+str(x+1)+'.png'
            
            fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(12,9))
            plt.subplots_adjust(left=.1, bottom=.05, \
                                right=.99, top=.95, wspace=.3, hspace=.4)
            for y in range(Ny):
                px, py = int(y/3), y%3
                axes[px,py].scatter(parmcmc[:Nmcmc,ipar,y,x], counter, s=.5)
                # for band in TABand:
                #     if band['label']==bname:
                #         if bex==0:
                #             axes[px,py].vlines(band['wave'],0,Nmcmc,colors='r')
                #         elif bex==1:
                #             axes[px,py].vlines(band['sigmaS'],0,Nmcmc,colors='r')
                #         elif bex==2:
                #             axes[px,py].vlines(band['sigmaL'],0,Nmcmc,colors='r')
                if fixed[ipar]=='F':
                    axes[px,py].vlines(chi2ini[ipar,y,x],0,Nmcmc,colors='r',label='Chi2')
                    if h5_sim.exists():
                        axes[px,py].vlines(simpar[ipar,y,x],0,Nmcmc,colors='g',label='True')
                
                axes[px,py].set_xlabel('Parameter value')
                axes[px,py].set_ylabel('MCMC counter')
                axes[px,py].set_title(
                    # fr'${parname[ipar]}({x+1},{y+1})\ [{meanpar[ipar,y,x]:.2f},{stdevpar[ipar,y,x]:.2f}]$')
                    fr'${parname[ipar]}({x+1},{y+1})={qpar[1,ipar,y,x]:.2f}$')
                axes[px,py].legend(loc='upper left')
            
            plt.savefig(filename)
        
        
            
        ## Parameter distribution (pix)
        ##------------------------------
        if Nmcmc<100:
            Nbin = Nmcmc
        elif Nmcmc<1000:
            Nbin = int(Nmcmc/10)
        else:
            Nbin = int(Nmcmc/100)
            
        for x in range(Nx):
            filename = path_fig+'pardist_'+m+'_x='+str(x+1)+'.png'
        
            fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(12,9))
            plt.subplots_adjust(left=.1, bottom=.05,
                                right=.99, top=.95, wspace=.3, hspace=.4)
            for y in range(Ny):
                px, py = int(y/3), y%3
                n, bins, patches = axes[px,py].hist(parmcmc[:Nmcmc,ipar,y,x],Nbin,
                                                    facecolor='g',alpha=.5,label='HB')
                if fixed[ipar]=='F':
                    axes[px,py].vlines(chi2ini[ipar,y,x],0,Nmcmc,colors='r',label='Chi2')
                    if h5_sim.exists():
                        axes[px,py].vlines(simpar[ipar,y,x],0,Nmcmc,colors='g',label='True')
                # for band in TABand:
                #     if band['label']=='Main 11.2':
                #         if bex==0:
                #             axes[px,py].vlines(band['wave'],0,Nmcmc,colors='r')
                #         elif bex==1:
                #             axes[px,py].vlines(band['sigmaS'],0,Nmcmc,colors='r')
                #         elif bex==2:
                #             axes[px,py].vlines(band['sigmaL'],0,Nmcmc,colors='r')
                axes[px,py].set_xlabel('Parameter value')
                axes[px,py].set_ylabel('Counts')
                axes[px,py].set_title(
                    # fr'${parname[ipar]}({x+1},{y+1})\ [{meanpar[ipar,y,x]:.2f},{stdevpar[ipar,y,x]:.2f}]$')
                    fr'${parname[ipar]}({x+1},{y+1})={qpar[1,ipar,y,x]:.2f}$')
                axes[px,py].legend(loc='upper left')
            
            plt.savefig(filename)
        
        
        
        ## Parameter track (par)
        ##-----------------------
        x,y = 1,6
        for j in range(int(Npar/9)+1):
            filename = path_fig+'partrack_'+m+'_ipar='+str(j*9+1)+'-'+str(j*9+9)+'.png'
            
            fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(12,9))
            plt.subplots_adjust(left=.1, bottom=.05, \
                                right=.99, top=.95, wspace=.3, hspace=.4)
        
            for jpar in range(9):
                jpar += j*9
                px, py = int(jpar%9/3), jpar%9%3
                if jpar<Npar:
                    axes[px,py].scatter(parmcmc[:Nmcmc,jpar,y,x], counter, s=.5)
        
                    ## Fitted (and after burning) param in magenta
                    axes[px,py].vlines(qpar[1,jpar,y,x],0,Nmcmc,colors='m',label='HB')
                    axes[px,py].vlines(qpar[0,jpar,y,x],0,Nmcmc,colors='m',linestyles='dashed')
                    axes[px,py].vlines(qpar[2,jpar,y,x],0,Nmcmc,colors='m',linestyles='dashed')
                    # axes[px,py].vlines(meanpar[jpar,y,x],0,Nmcmc,colors='m')
                    # axes[px,py].vlines(meanpar[jpar,y,x]-stdevpar[jpar,y,x],
                    #                    0,Nmcmc,colors='m',linestyles='dashed')
                    # axes[px,py].vlines(meanpar[jpar,y,x]+stdevpar[jpar,y,x],
                    #                    0,Nmcmc,colors='m',linestyles='dashed')
        
                    ## Initial param in red
                    if fixed[jpar]=='F':
                        # axes[px,py].vlines(parini[jpar],0,Nmcmc,colors='r')
                        axes[px,py].vlines(chi2ini[jpar,y,x],0,Nmcmc,colors='r',label='Chi2')
                        if h5_sim.exists():
                            axes[px,py].vlines(simpar[jpar,y,x],0,Nmcmc,colors='g',label='True')
                    
                    axes[px,py].set_xlabel('Parameter value')
                    axes[px,py].set_ylabel('MCMC counter')
                    axes[px,py].set_title(
                        # fr'${parname[jpar]}({x+1},{y+1})\ [{meanpar[jpar,y,x]:.2f},{stdevpar[jpar,y,x]:.2f}]$')
                        fr'${parname[jpar]}({x+1},{y+1})={qpar[1,jpar,y,x]:.2f}$')
                    axes[px,py].legend(loc='upper left')
            
            plt.savefig(filename)
        
        
        
        ## ACFs
        ##------
        ihyp = 23 # lnRband5 6.2 (1)
        icorr = 1867
        hyparname = read_hdf5(filout, 'Hyperparameter label')
        
        autocor_mu = read_hdf5(filout, 'Autocorrelation function for mu('+hyparname[ihyp]+')')
        autocor_sig = read_hdf5(filout, 'Autocorrelation function for sig('+hyparname[ihyp]+')')
        mcmc_hyp = np.arange(len(autocor_mu))
        
        autocor_corr = read_hdf5(filout, 'Autocorrelation function for corr('+str(icorr)+')')
        # autocor_corr = read_hdf5(filout, 'Autocorrelation function for corr('+str(icorr)+') after burn-in')
        mcmc_corr = np.arange(len(autocor_corr))
        
        fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(8,8))
        plt.subplots_adjust(left=.1, bottom=.05, \
                            right=.99, top=.95, wspace=.3, hspace=.4)
        
        axes[0].plot(mcmc_hyp, autocor_mu)
        axes[0].set_xlabel('Lag')
        axes[0].set_ylabel('ACF of mu('+hyparname[ihyp]+')')
        axes[1].plot(mcmc_hyp, autocor_sig)
        axes[1].set_xlabel('Lag')
        axes[1].set_ylabel('ACF of sig('+hyparname[ihyp]+')')
        axes[2].plot(mcmc_corr, autocor_corr)
        axes[2].set_xlabel('Lag')
        axes[2].set_ylabel('ACF of corr'+str(icorr))
        xmin, xmax = 0, 50
        axes[0].set_xlim((xmin,xmax))
        axes[1].set_xlim((xmin,xmax))
        axes[2].set_xlim((xmin,xmax))
        
        plt.savefig(path_fig+'partrack_'+m+'_ACFs.png')

# plt.show()

print('>>> Coucou show_par [Done] <<<')
