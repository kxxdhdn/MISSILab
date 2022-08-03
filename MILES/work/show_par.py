#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Tracking parameter distribution during Gibbs sampling

"""

import os, pathlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

## rapyuta
from rapyuta.maths import ij2icorr
from rapyuta.inout import read_hdf5, h5ext
from rapyuta.plots import pplot, plotool

## local
from auxil import croot, plotname

## Path
##------
path_out = croot+'/../out/'
path_fig = path_out+'Figures/'
os.makedirs(path_fig, exist_ok=True)
filobs = path_out+'observation_MIR'
spec_name = read_hdf5(filobs, 'spectrum labels')
filsim = path_out+'simulation_MIR'
h5_sim = pathlib.Path(filsim+h5ext)
if h5_sim.exists():
    simpar = read_hdf5(filsim, 'Simulated parameter value') # True par
    simcal = read_hdf5(filsim, 'Simulated ln(1+delta)') # True ln1pd
    # simcal = read_hdf5(path_out+'presimu_bb', 'Quantiles of ln(1+delta)')[1,:,0,0] # True ln1pd
chi2_out = pathlib.Path(path_out+'fit_chi2'+h5ext)
if chi2_out.exists():
    chi2ini = read_hdf5(path_out+'fit_chi2', 'Best fitted parameter value') # Chi2 result
h5_analysis = path_out+'input_analysis'

## Read h5 file
mode = ['bb','hb']
calib = True

##---------
## Indices
##---------
ipar = 6 # lnT2
## y = 0 (c1_SN10), 1 (c10_SN10), 2 (c1_SN100), 3 (c10_SN100)
y = 1
x = 0

## Multivariate plot
impar = [1, 2, 4, 6, 78] # 78
impar2 = [1, 2, 4, 6, 78]
suffix = 'cont' # [Lstar, Fcont, Ibandref] - [Lstar - Fcont, Ibandref]
# impar = [1, 2, 4, 6, 78]
# impar2 = [0]
# suffix = 'avext' # [Lstar, Fcont, Ibandref] - [Av]
# impar = [11, 38, 58, 94]
# impar2 = [0]
# suffix = 'avrat' # [Rline, Rband] - [Av]
# impar = [11, 14, 58, 62, 66]
# impar2 = [11, 14, 58, 62, 66]
# suffix = 'rat' # [Rline, Rband] - [Rline, Rband]
Nmult = len(impar) # X - width
Nmult2 = len(impar2) # Y - height
multivar = input('Multivariate plot (y/n)? ')

for m in mode:
    if m=='chi2':
        MODE = 'Chi2'
    elif m=='bb':
        MODE = 'non-HB'
        ibm = 0
    elif m=='hb':
        MODE = 'HB'
        ibm = 1
    filout = path_out+'fit_'+m
    filmcmc = path_out+'parlog_fit_'+m

    h5_out = pathlib.Path(filout+h5ext)
    if h5_out.exists():
            
        # Nmcmc = read_hdf5(filmcmc, 'Length of MCMC')[0]
        Nmcmc = read_hdf5(filmcmc, 'Last index')[0]
        counter = np.arange(Nmcmc)
        # t_end = Nmcmc
        # t_burnin = int(t_end*.3) - 1
        t_end = read_hdf5(h5_analysis, 't_end')[ibm]
        t_burnin = read_hdf5(h5_analysis, 't_burnin')[ibm]
        parname = read_hdf5(filmcmc, 'Parameter label')
        parmcmc = read_hdf5(filmcmc, 'Parameter values')
        instr = read_hdf5(filobs, 'spectroscopic module labels')
        
        Npar, Ny, Nx = parmcmc.shape[1:4]
        meanpar = read_hdf5(filout, 'Mean of parameter value')
        # mu = np.mean(parmcmc[t_burnin:t_end,:,:,:], axis=0)
        # print(mu-meanpar)
        stdevpar = read_hdf5(filout, 'Sigma of parameter value')
        # sig = np.std(parmcmc[t_burnin:t_end,:,:,:], axis=0)
        # print(sig-stdevpar)
        # exit()
        qpar = read_hdf5(filout, 'Quantiles of parameter value')

        if calib:
            ln1pdmcmc = read_hdf5(filmcmc, 'ln(1+delta)')
            Ncal = ln1pdmcmc.shape[1]
            meanln1pd = read_hdf5(filout, 'Mean of ln(1+delta)')
            stdevln1pd = read_hdf5(filout, 'Sigma of ln(1+delta)')
            qln1pd = read_hdf5(filout, 'Quantiles of ln(1+delta)')
        
        filmod = path_out+'input_model'
        fixed = read_hdf5(filmod, 'parinfo fixed')
        parini = read_hdf5(filmod, 'parinfo value')
        labB = read_hdf5(filmod, 'label band')

        if multivar=='y':
            
            ## Multivariate plot
            ##-------------------
            nampar = []
            width_ratios = [1]
            for i in range(Nmult):
                nampar.append(parname[impar[i]])
                width_ratios.append(2)
            nampar2 = []
            height_ratios = []
            for i in range(Nmult2):
                nampar2.append(parname[impar2[i]])
                height_ratios.append(2)
            height_ratios.append(1)
            print(spec_name[y,x], MODE)
            print('Multivariate : '+str(nampar)+' - '+str(nampar2))
            space_width = .6
            fxsize = (Nmult+.4) * (space_width+2)
            fysize = (Nmult2+.4) * (space_width+2)
            # print(fxsize, fysize)
            
            p = plotool(Nmult2+1, Nmult+1, figsize=(fxsize,fysize),
                        # sharex=True, sharey=True,
                        gridspec_kw={ 'height_ratios':height_ratios,
                                      'width_ratios':width_ratios })
            p.set_fig(hspace=space_width, wspace=space_width)
            
            ## blank (lower left) frame
            ##--------------------------
            p.set_ax([Nmult2,0], xlabel='Counts', ylabel='Counts',
                     left=False,bottom=False,labelleft=False,labelbottom=False,
                     xsize=15,ysize=15,xtksize=15,ytksize=15,)
            for spine in p.ax.spines.values():
                spine.set_edgecolor('None')
            # p.fig.delaxes(p.axes[Nmult2,0])
            
            for px in range(Nmult2):
                for py in range(Nmult):
            
                    xmin = min(parmcmc[:Nmcmc,impar[py],y,x])
                    xmax = max(parmcmc[:Nmcmc,impar[py],y,x])
                    if h5_sim.exists():
                        xmin = min( np.append(xmin,simpar[impar[py],y,x]) )
                        xmax = max( np.append(xmax,simpar[impar[py],y,x]) )
                    if chi2_out.exists():
                        xmin = min( np.append(xmin,chi2ini[impar[py],y,x]) )
                        xmax = max( np.append(xmax,chi2ini[impar[py],y,x]) )
                    xgap = xmax - xmin
                    xmin += -xgap* .05
                    xmax += xgap* .05
                    
                    ymin = min(parmcmc[:Nmcmc,impar2[px],y,x])
                    ymax = max(parmcmc[:Nmcmc,impar2[px],y,x])
                    if h5_sim.exists():
                        ymin = min( np.append(ymin,simpar[impar2[px],y,x]) )
                        ymax = max( np.append(ymax,simpar[impar2[px],y,x]) )
                    if chi2_out.exists():
                        ymin = min( np.append(ymin,chi2ini[impar2[px],y,x]) )
                        ymax = max( np.append(ymax,chi2ini[impar2[px],y,x]) )
                    ygap = ymax - ymin
                    ymin += -ygap* .05
                    ymax += ygap* .05
                    
                    ## histogram (along X)
                    ##---------------------
                    xlab = plotname(nampar[py])
                    # if py==0:
                    #     ylab = 'Counts'
                    # else:
                    ylab = None
                        
                    p.set_ax([Nmult2,py+1],xlabel=xlab,ylabel=ylab,
                             xlim=(xmin,xmax), labelrotation=45,
                             xsize=15,ysize=15,xtksize=15,ytksize=15,)
                    count, bins, patches = p.ax.hist(
                        parmcmc[t_burnin:t_end,impar[py],y,x], bins='sqrt', log=False,
                        facecolor='c', edgecolor='None', alpha=.5,)
                    ## Simu param in green
                    if h5_sim.exists():
                        p.ax.axvline(simpar[impar[py],y,x],c='g',lw=4,alpha=.5,zorder=100)
                    ## Chi2 param in red
                    if chi2_out.exists():
                        p.ax.axvline(chi2ini[impar[py],y,x],c='r',zorder=100)
                    ## Fitted (and after burn-in) param in magenta
                    p.ax.axvline(qpar[1,impar[py],y,x],c='m',zorder=100)
                    p.ax.axvline(qpar[0,impar[py],y,x],c='m',ls='dashed',zorder=100)
                    p.ax.axvline(qpar[2,impar[py],y,x],c='m',ls='dashed',zorder=100)
                    # p.append_handles()
                    # p.set_legend(loc='extra upper left',locext=.4,fontsize=15,framealpha=0)
                    # p.reset_handles()
            
                    ## histogram (along Y)
                    ##---------------------
                    ylab = plotname(nampar2[px])
                    # if px==Nmult2-1:
                    #     xlab = 'Counts'
                    # else:
                    xlab = None
                        
                    p.set_ax([px,0],xlabel=xlab,ylabel=ylab,
                             ylim=(ymin,ymax), labelrotation=45,
                             xsize=15,ysize=15,xtksize=15,ytksize=15,)
                    count, bins, patches = p.ax.hist(
                        parmcmc[t_burnin:t_end,impar2[px],y,x], bins='sqrt', log=False,
                        facecolor='c', edgecolor='None', alpha=.5,
                        orientation="horizontal")
                    ## Simu param in green
                    if h5_sim.exists():
                        p.ax.axhline(simpar[impar2[px],y,x],c='g',lw=4,alpha=.5,zorder=100)
                    ## Chi2 param in red
                    if chi2_out.exists():
                        p.ax.axhline(chi2ini[impar2[px],y,x],c='r',zorder=100)
                    ## Fitted (and after burn-in) param in magenta
                    p.ax.axhline(qpar[1,impar2[px],y,x],c='m',zorder=100)
                    p.ax.axhline(qpar[0,impar2[px],y,x],c='m',ls='dashed',zorder=100)
                    p.ax.axhline(qpar[2,impar2[px],y,x],c='m',ls='dashed',zorder=100)
            
                    ## cloud scattering
                    ##------------------
                    p.set_ax([px,py+1], labelrotation=45,
                             labelleft=False,labelbottom=False,
                             xlim=(xmin,xmax), ylim=(ymin,ymax),
                             xsize=15,ysize=15,xtksize=15,ytksize=15,)
                    p.plot(parmcmc[t_burnin:t_end,impar[py],y,x],
                           parmcmc[t_burnin:t_end,impar2[px],y,x],
                           fmt='s',c='c',marker='.',markersize=10,
                           markeredgecolor='None',alpha=0.2)
                    ## Simu param in green
                    if h5_sim.exists():
                        p.plot(simpar[impar[py],y,x], simpar[impar2[px],y,x],
                               fmt='s',c='g',marker='*',markersize=20,zorder=100)
                    ## Chi2 param in red
                    if chi2_out.exists():
                        p.plot(chi2ini[impar[py],y,x], chi2ini[impar2[px],y,x],
                               fmt='s',c='r',marker='*',markersize=10,zorder=100)
                    
                    ## lower triangle for same parameters
                    if impar[py] in impar2:
                        if impar[py]>=impar2[px]:
                            if py+1>=px:
                                p.fig.delaxes(p.axes[px,py+1])
                            else:
                                p.fig.delaxes(p.axes[py+1,px])
                
            filename = path_fig+'multivariate_'+m+'_'+spec_name[y,x]+'_'+suffix+'.png'
            p.save(filename, transparent=True, figtight=True)

        else:
            
            '''
            ## Parameter track (pix)
            ##-----------------------
            print(parname[ipar], MODE)
            
            for x in range(Nx):
                for y in range(Ny):
            
                    parlab = plotname(parname[ipar])

                    p = plotool(2, 1, figsize=(9,8), gridspec_kw={'height_ratios':[3,1]})
                    p.set_fig(hspace=.1)
            
                    ## walk
                    p.set_ax([0,0], xlabel='Parameter value', ylabel='MCMC counter',
                             xtktoggle=True, xsize=15, ysize=15, xtksize=15, ytksize=15,)
                    p.plot(parmcmc[:Nmcmc,ipar,y,x], counter,
                           fmt='s', c='c', marker='.', markersize=5, label=parlab)
                    ## Simu param in green
                    if h5_sim.exists():
                        p.ax.axvline(simpar[ipar,y,x],c='g',lw=4,alpha=.5,
                                     label=fr'True: ${simpar[ipar,y,x]:.3f}$',zorder=100)
                    ## Chi2 param in red
                    if chi2_out.exists():
                        p.ax.axvline(chi2ini[ipar,y,x],c='r',
                                     label=fr'Chi2: ${chi2ini[ipar,y,x]:.3f}$',zorder=100)
                    ## Fitted (and after burn-in) param in magenta
                    p.ax.axvline(qpar[1,ipar,y,x],c='m',
                                 label=MODE+fr': ${qpar[1,ipar,y,x]:.3f}$',zorder=100)
                    p.ax.axvline(qpar[0,ipar,y,x],c='m',ls='dashed',zorder=100)
                    p.ax.axvline(qpar[2,ipar,y,x],c='m',ls='dashed',zorder=100)
                    # p.ax.axvline(meanpar[ipar,y,x],c='m',zorder=100)
                    # p.ax.axvline(meanpar[ipar,y,x]-stdevpar[ipar,y,x],
                    #            c='m',ls='dashed',zorder=100)
                    # p.ax.axvline(meanpar[ipar,y,x]+stdevpar[ipar,y,x],
                    #            c='m',ls='dashed',zorder=100)
                    p.append_handles()
                    p.set_legend(loc='extra upper left', locext=.4, fontsize=15, framealpha=0)
                    p.reset_handles()
                    p.ax.invert_yaxis()
            
                    ## dist (burn-in)
                    p.set_ax([1,0], xlabel='Parameter value', ylabel='MCMC counts',
                             xsize=15, ysize=15, xtksize=15, ytksize=15,)
                    count, bins, patches = p.ax.hist(
                        parmcmc[t_burnin:t_end,ipar,y,x], bins='sqrt', log=False,
                        facecolor='c', edgecolor='None', alpha=.5, label='(burn-in)')
                    # ## Simu param in green
                    # if h5_sim.exists():
                    #     p.ax.axvline(simpar[ipar,y,x],c='g',lw=4,alpha=.5,
                    #                  label='True',zorder=100)
                    # ## Chi2 param in red
                    # if chi2_out.exists():
                    #     p.ax.axvline(chi2ini[ipar,y,x],c='r',label='Chi2',zorder=100)
                    # ## Fitted (and after burn-in) param in magenta
                    # p.ax.axvline(qpar[1,ipar,y,x],c='m',label=MODE,zorder=100)
                    # p.ax.axvline(qpar[0,ipar,y,x],c='m',ls='dashed',zorder=100)
                    # p.ax.axvline(qpar[2,ipar,y,x],c='m',ls='dashed',zorder=100)
                    p.append_handles()
                    p.set_legend(loc='extra upper left', locext=.4, fontsize=15, framealpha=0)
            
                    
                    filename = path_fig+'par_'+m+'_'+parname[ipar]+'_'+spec_name[y,x]+'.png'
                    p.save(filename, transparent=True, figtight=True)
            '''
            
            ## Parameter track (par)
            ##-----------------------
            print(spec_name[y,x], MODE)
            
            jj = 0
            for j in range(Npar):
            
                parlab = plotname(parname[j])
                
                if fixed[j]=='F':
                    jj += 1
            
                    p = plotool(2, 1, figsize=(9,8), gridspec_kw={'height_ratios':[3,1]})
                    p.set_fig(hspace=.1)
            
                    xmin = min(parmcmc[:Nmcmc,j,y,x])
                    xmax = max(parmcmc[:Nmcmc,j,y,x])
                    if h5_sim.exists():
                        xmin = min( np.append(xmin,simpar[j,y,x]) )
                        xmax = max( np.append(xmax,simpar[j,y,x]) )
                    if chi2_out.exists():
                        xmin = min( np.append(xmin,chi2ini[j,y,x]) )
                        xmax = max( np.append(xmax,chi2ini[j,y,x]) )
                    xgap = xmax - xmin
                    xmin += -xgap* .05
                    xmax += xgap* .05
            
                    ## walk
                    p.set_ax([0,0], xlabel='Parameter value', ylabel='MCMC counter',
                             xlim=(xmin,xmax),
                             xtktoggle=True, xsize=15, ysize=15, xtksize=15, ytksize=15,)
                    p.plot(parmcmc[:Nmcmc,j,y,x], counter,
                           fmt='s', c='c', marker='.', markersize=5, label=parlab)
                    ## Simu param in green
                    if h5_sim.exists():
                        p.ax.axvline(simpar[j,y,x],c='g',lw=4,alpha=.5,
                                     label=fr'True: ${simpar[j,y,x]:.3f}$',zorder=100)
                    ## Chi2 param in red
                    if chi2_out.exists():
                        p.ax.axvline(chi2ini[j,y,x],c='r',
                                     label=fr'Chi2: ${chi2ini[j,y,x]:.3f}$',zorder=100)
                    ## Fitted (and after burn-in) param in magenta
                    p.ax.axvline(qpar[1,j,y,x],c='m',
                                 label=MODE+fr': ${qpar[1,j,y,x]:.3f}$',zorder=100)
                    p.ax.axvline(qpar[0,j,y,x],c='m',ls='dashed',zorder=100)
                    p.ax.axvline(qpar[2,j,y,x],c='m',ls='dashed',zorder=100)
                    # p.ax.axvline(meanpar[j,y,x],c='m',zorder=100)
                    # p.ax.axvline(meanpar[j,y,x]-stdevpar[j,y,x],c='m',ls='dashed',zorder=100)
                    # p.ax.axvline(meanpar[j,y,x]+stdevpar[j,y,x],c='m',ls='dashed',zorder=100)
                    p.append_handles()
                    p.set_legend(loc='extra upper left', locext=.4, fontsize=15, framealpha=0)
                    p.reset_handles()
                    # p.ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
                    p.ax.invert_yaxis()
            
                    ## dist (burn-in)
                    p.set_ax([1,0],xlabel='Parameter value',ylabel='MCMC counts',
                             xlim=(xmin,xmax),
                             xsize=15,ysize=15,xtksize=15,ytksize=15,)
                    count, bins, patches = p.ax.hist(
                        parmcmc[t_burnin:t_end,j,y,x], bins='sqrt', log=False,
                        facecolor='c', edgecolor='None', alpha=.5, label='After burn-in')
                    ## Simu param in green
                    if h5_sim.exists():
                        p.ax.axvline(simpar[j,y,x],c='g',lw=4,alpha=.5,zorder=100)
                    ## Chi2 param in red
                    if chi2_out.exists():
                        p.ax.axvline(chi2ini[j,y,x],c='r',zorder=100)
                    ## Fitted (and after burn-in) param in magenta
                    p.ax.axvline(qpar[1,j,y,x],c='m',zorder=100)
                    p.ax.axvline(qpar[0,j,y,x],c='m',ls='dashed',zorder=100)
                    p.ax.axvline(qpar[2,j,y,x],c='m',ls='dashed',zorder=100)
                    p.append_handles()
                    p.set_legend(loc='extra upper left',locext=.4,fontsize=15,framealpha=0)
                
                    filename = path_fig+'par_'+m+'_'+spec_name[y,x]+'_'+str(jj)+'_'+str(j+1)+'.png'
                    p.save(filename, transparent=True, figtight=True)
            
            
            ## Parameter track (calib)
            ##-------------------------
            if calib:
                for j in range(1,Ncal):
                    print(instr[j])
            
                    p = plotool(2, 1, figsize=(8,8), gridspec_kw={'height_ratios':[3,1]})
                    p.set_fig(hspace=.1)
            
                    xmin = min(np.exp(ln1pdmcmc[:Nmcmc,j,y,x]))
                    xmax = max(np.exp(ln1pdmcmc[:Nmcmc,j,y,x]))
                    if h5_sim.exists():
                        xmin = min( np.append(xmin,np.exp(simcal[j])) )
                        xmax = max( np.append(xmax,np.exp(simcal[j])) )
                    xgap = xmax - xmin
                    xmin += -xgap* .05
                    xmax += xgap* .05
                
                    ## walk
                    p.set_ax([0,0],xlabel='Calibration factor',ylabel='MCMC counter',
                             xlim=(xmin,xmax),
                             xtktoggle=True,xsize=15,ysize=15,xtksize=15,ytksize=15,)
                    p.plot(np.exp(ln1pdmcmc[:Nmcmc,j,y,x]),counter,
                           fmt='s',c='c',marker='.',markersize=5,label=instr[j],)
                    ## Simu param in green
                    if h5_sim.exists():
                        p.ax.axvline(np.exp(simcal[j]),c='g',lw=4,alpha=.5,
                                     label=fr'True: ${np.exp(simcal[j]):.3f}$',zorder=100)
                    ## Fitted (and after burn-in) param in magenta
                    p.ax.axvline(np.exp(qln1pd[1,j,y,x]),c='m',
                                 label=MODE+fr': ${np.exp(qln1pd[1,j,y,x]):.3f}$',zorder=100)
                    p.ax.axvline(np.exp(qln1pd[0,j,y,x]),c='m',ls='dashed',zorder=100)
                    p.ax.axvline(np.exp(qln1pd[2,j,y,x]),c='m',ls='dashed',zorder=100)
                    # p.ax.axvline(np.exp(meanln1pd[j,y,x]),c='m',label=MODE,zorder=100)
                    # p.ax.axvline(np.exp(meanln1pd[j,y,x]-stdevln1pd[j,y,x]),c='m',ls='dashed',zorder=100)
                    # p.ax.axvline(np.exp(meanln1pd[j,y,x]+stdevln1pd[j,y,x]),c='m',ls='dashed',zorder=100)
                    p.append_handles()
                    p.set_legend(loc='extra upper left',locext=.4,fontsize=15,framealpha=0)
                    p.reset_handles()
                    p.ax.invert_yaxis()
                
                    ## dist (burn-in)
                    p.set_ax([1,0],xlabel='Calibration factor',ylabel='MCMC counts',
                             xlim=(xmin,xmax),
                             xsize=15,ysize=15,xtksize=15,ytksize=15,)
                    count, bins, patches = p.ax.hist(
                        np.exp(ln1pdmcmc[t_burnin:t_end,j,y,x]),bins='sqrt',log=False,
                        facecolor='c',edgecolor='None',alpha=.5,label='After burn-in')
                    ## Simu param in green
                    if h5_sim.exists():
                        p.ax.axvline(np.exp(simcal[j]),c='g',lw=4,alpha=.5,zorder=100)
                    ## Fitted (and after burn-in) param in magenta
                    p.ax.axvline(np.exp(qln1pd[1,j,y,x]),c='m',zorder=100)
                    p.ax.axvline(np.exp(qln1pd[0,j,y,x]),c='m',ls='dashed',zorder=100)
                    p.ax.axvline(np.exp(qln1pd[2,j,y,x]),c='m',ls='dashed',zorder=100)
                    p.append_handles()
                    p.set_legend(loc='extra upper left',locext=.4,fontsize=15,framealpha=0)
                
                    filename = path_fig+'calib_'+m+'_'+spec_name[y,x]+'_'+instr[j]+'.png'
                    p.save(filename, transparent=True, figtight=True)
        

print('>>> Coucou show_par [Done] <<<')
