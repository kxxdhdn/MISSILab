#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Tracking parameter distribution during Gibbs sampling

"""

import os, pathlib
import numpy as np
from scipy.stats import kde
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

## rapyuta
from rapyuta.inout import read_hdf5, h5ext
from rapyuta.plots import pplot, plotool

## local
from auxil import croot, plotcorr, calcorr

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
    chi2err = read_hdf5(path_out+'fit_chi2', 'Best fitted parameter error')
h5_analysis = path_out+'input_analysis'

## Read h5 file
mode = ['bb','hb']
calib = True

## Prepare correlations
corrname = [
    'I3.3_ov_I11.3',
    'I3.4_ov_I11.3',
    'I6.2_ov_I11.3',
    'I7.7_ov_I11.3',
    'I12.7_ov_I11.3',
]
## 11.3, ...
indpar = [82, 38, 42,46, 50,54, 58,62,66, 86,90]

xcorr = [1,0,2,3]
ycorr = [1,0,2,3,4]
suffix = 'c0'

# xcorr = [0]
# ycorr = [3]
# suffix = 'c14'

# xcorr = [1]
# ycorr = [2]
# suffix = 'c23'

# xcorr = [1]
# ycorr = [4]
# suffix = 'c25'

# xcorr = [2]
# ycorr = [3]
# suffix = 'c34'

Nmult = len(xcorr) # X - width
Nmult2 = len(ycorr) # Y - height

den = 0 # density plot / scatter plot

group = 4

if group==1:
    i1,i2 = 0,1
    gsuf = '_g1'
elif group==2:
    i1,i2 = 1,2
    gsuf = '_g2'
elif group==3:
    i1,i2 = 2,3
    gsuf = '_g3'
elif group==4:
    i1,i2 = 3,4
    gsuf = '_g4'
else:
    i1,i2 = 0,4
    gsuf = ''

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

        ## Prepare data
        ##--------------
        colors = np.empty((Ny,Nx), dtype=('<U30'))
        zorders = np.zeros((Ny,Nx))
        for y in range(Ny):
            for x in range(Nx):
                if y==0:
                    colors[y,x] = 'c'
                    zorders[y,x] = 500
                elif y==1:
                    colors[y,x] = 'g'
                    zorders[y,x] = 0
                elif y==2:
                    colors[y,x] = 'm'
                    zorders[y,x] = 500
                elif y==3:
                    colors[y,x] = 'y'
                    zorders[y,x] = 500
        colorall = colors.reshape(Ny*Nx)
        zall = zorders.reshape(Ny*Nx)

        simcorr = []
        chi2corr = []
        chi2correrr = []
        corrmcmc = []
        for cnam in corrname:
            simcorr0 = np.empty(Ny*Nx).reshape(Ny,Nx)
            chi2corr0 = np.empty(Ny*Nx).reshape(Ny,Nx)
            chi2correrr0 = np.empty(Ny*Nx).reshape(Ny,Nx)
            corrmcmc0 = np.empty(Nmcmc*Ny*Nx).reshape(Nmcmc,Ny,Nx)
            for y in range(Ny):
                for x in range(Nx):
                    for i in range(Nmcmc):
                        corrmcmc0[i,y,x] = calcorr(parmcmc[i,:,y,x], cnam, indpar)[0]
                    simcorr0[y,x] = calcorr(simpar[:,y,x], cnam, indpar)[0]
                    chi2corr0[y,x], chi2correrr0[y,x] = calcorr(
                        chi2ini[:,y,x], cnam, indpar, chi2err[:,y,x], chi2=True)
            simcorr.append( simcorr0.reshape(Ny*Nx) )
            chi2corr.append( chi2corr0.reshape(Ny*Nx) )
            chi2correrr.append( chi2correrr0.reshape(Ny*Nx) )
            corrmcmc.append( corrmcmc0.reshape(Nmcmc,Ny*Nx) )
        simcorr = np.array(simcorr)
        chi2corr = np.array(chi2corr)
        chi2correrr = np.array(chi2correrr)
        corrmcmc = np.array(corrmcmc)
        
        ## Correlation plot
        ##------------------
        namcorr = []
        width_ratios = []
        for i in range(Nmult):
            namcorr.append(corrname[xcorr[i]])
            width_ratios.append(2)
        namcorr2 = []
        height_ratios = []
        for i in range(Nmult2):
            namcorr2.append(corrname[ycorr[i]])
            height_ratios.append(2)
        # height_ratios.append(1)
        print(MODE)
        print('Multivariate : '+str(namcorr)+' - '+str(namcorr2))
        space_width = .1
        fxsize = (Nmult+.4) * (space_width+2)
        fysize = (Nmult2+.4) * (space_width+2)
        # print(fxsize, fysize)
        
        p = plotool(Nmult2, Nmult, figsize=(fxsize,fysize),)
                    # sharex=True, sharey=True,
                    # gridspec_kw={ 'height_ratios':height_ratios,
                    #               'width_ratios':width_ratios })
        p.set_fig(hspace=space_width, wspace=space_width)
        
        for px in range(Nmult2):
            for py in range(Nmult):

                ## Prepare data
                ##--------------
                
                ## Prepare plot
                ##--------------
                # xmin,xmax,ymin,ymax=None,None,None,None
                xmin = min(corrmcmc[xcorr[py],t_burnin:t_end,i1*Nx:i2*Nx].flatten())
                xmax = max(corrmcmc[xcorr[py],t_burnin:t_end,i1*Nx:i2*Nx].flatten())
                if h5_sim.exists():
                    xmin = min( np.append(xmin,simcorr[xcorr[py],i1*Nx:i2*Nx]) )
                    xmax = max( np.append(xmax,simcorr[xcorr[py],i1*Nx:i2*Nx]) )
                # if chi2_out.exists():
                #     xmin = min( np.append(xmin,chi2corr[xcorr[py],i1*Nx:i2*Nx]) )
                #     xmax = max( np.append(xmax,chi2corr[xcorr[py],i1*Nx:i2*Nx]) )
                # xgap = xmax - xmin
                # xmin += -xgap* .05
                # xmax += xgap* .05
                xmin *= 0.8
                xmax *= 1.2
                
                ymin = min(corrmcmc[ycorr[px],t_burnin:t_end,i1*Nx:i2*Nx].flatten())
                ymax = max(corrmcmc[ycorr[px],t_burnin:t_end,i1*Nx:i2*Nx].flatten())
                if h5_sim.exists():
                    ymin = min( np.append(ymin,simcorr[ycorr[px],i1*Nx:i2*Nx]) )
                    ymax = max( np.append(ymax,simcorr[ycorr[px],i1*Nx:i2*Nx]) )
                # if chi2_out.exists():
                #     ymin = min( np.append(ymin,chi2corr[ycorr[px],i1*Nx:i2*Nx]) )
                #     ymax = max( np.append(ymax,chi2corr[ycorr[px],i1*Nx:i2*Nx]) )
                # ygap = ymax - ymin
                # ymin += -ygap* .05
                # ymax += ygap* .05
                ymin *= 0.8
                ymax *= 1.2
                
                # xmin,ymin = 5e-2, 5e-2
                # xmax,ymax = 2e1, 2e1
                # xmin,ymin = None, None
                # xmax,ymax = None, None

                xlab, ylab = None, None
                labelleft, labelbottom = False, False
                if px==Nmult2-1:
                    xlab = plotcorr(namcorr[py])
                    labelbottom = True
                if py==0:
                    ylab = plotcorr(namcorr2[px])
                    labelleft = True

                ## cloud scattering
                ##------------------
                p.set_ax([px,py], labelrotation=45,
                         xlabel=xlab,ylabel=ylab,
                         labelleft=labelleft,labelbottom=labelbottom,
                         xlim=(xmin,xmax), ylim=(ymin,ymax),
                         xlog=1,ylog=1,xtkform='log_sci',ytkform='log_sci',
                         xsize=15,ysize=15,xtksize=15,ytksize=15,)
                # for i in range(Ny*Nx):
                for i in range(i1*Nx,i2*Nx):
                    x = corrmcmc[xcorr[py],t_burnin:t_end,i]
                    y = corrmcmc[ycorr[px],t_burnin:t_end,i]
                    if den:
                        ## Density plot
                        nbins = 20
                        k = kde.gaussian_kde([x,y])
                        xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
                        zi = k(np.vstack([xi.flatten(), yi.flatten()]))
                        p.ax.pcolormesh(xi, yi, zi.reshape(xi.shape),
                                        shading='gouraud', cmap=plt.cm.BuPu)
                    else:
                        ## Scatter plot
                        p.plot(x,y,
                               fmt='s',c=colorall[i],marker='.',markersize=5,
                               markeredgecolor='None',alpha=0.1,zorder=zall[i])
                    ## Simu param in black
                    if h5_sim.exists():
                        p.plot(simcorr[xcorr[py],i], simcorr[ycorr[px],i],
                               fmt='s',c='k',marker='.',markersize=10,
                               markeredgecolor='None',alpha=1,zorder=1000)
                    ## Chi2 param in red
                    if chi2_out.exists():
                        p.plot(chi2corr[xcorr[py],i], chi2corr[ycorr[px],i],
                               yerr=chi2correrr[ycorr[px],i], xerr=chi2correrr[xcorr[py],i],
                               fmt='s',c='grey',ec='grey',marker='.',markersize=5,# capsize=2,
                               markeredgecolor='None',alpha=1,zorder=1001)

                ## lower triangle for same parameters
                if xcorr[py] in ycorr:
                    if xcorr[py]>=ycorr[px]:
                        if py>=px:
                            p.fig.delaxes(p.axes[px,py])
                        else:
                            p.fig.delaxes(p.axes[py,px])
            
        filename = path_fig+'corr_'+suffix+'_'+m+gsuf+'.png'
        p.save(filename, transparent=True, figtight=True)
        

print('>>> Coucou correlations [Done] <<<')
