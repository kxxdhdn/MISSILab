#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This is the visualization of correlations

"""

import os, pathlib
import numpy as np
from astropy.io import ascii
from scipy.optimize import curve_fit
from matplotlib.ticker import ScalarFormatter, NullFormatter
import matplotlib.colors as mcolors

## rapyuta
from rapyuta.arrays import closest
from rapyuta.inout import read_hdf5
from rapyuta.plots import pplot

## local
from auxil import croot
from lib_corr import (clist, mlist, labels, fnames,
                      values_chi2, errors_chi2, limits_chi2,
                      values_bb, errors_bb, limits_bb,
                      values_hb, errors_hb, limits_hb)


def func(x, a, b):
    '''
    f(x) = a * x + b
    '''
    return a * x + b

def non_corr_df1(x, xerr):
    '''
    f1(x) = exp(x)
    df1 = exp(x) * dx
    '''
    return xerr * np.exp(x)

def non_corr_df2(x, y, xerr, yerr):
    '''
    f2(x,y) = f1(x)*f1(y) = exp(x+y)
    df2 = exp(x+y) * sqrt(dx^2 + dy^2) (non-correlated)
    '''
    return np.sqrt(xerr**2 + yerr**2) * np.exp(x+y)


mode = ['chi2', 'bb', 'hb']


## Path
##------
path_out = croot+'../out/'
path_fig = path_out+'Figures/'
if not os.path.exists(path_fig):
    os.makedirs(path_fig)
h5_obs = path_out+'observation_MIR'
h5_model = path_out+'input_model'
h5_analysis = path_out+'input_analysis'

## Read csv
##----------
csv_model = path_out+'model_PAH'
ds = ascii.read(csv_model+'.csv')
I33_I113 = ds['col1']
I77_I113 = ds['col2']
fion = ds['col3']
Ncmin = ds['col4']
hnu = ds['col5']
G0 = ds['col6']
G0_nT05 = ds['col7']
ISRFtype = ds['col8']
Ngrid = len(I33_I113)

## Read fit
##----------
spec_name = read_hdf5(h5_obs, 'spectrum labels')

for m in mode:
    filout = path_out+'fit_'+m
    h5_out = pathlib.Path(filout+'.h5')
    if h5_out.exists():
        print(m)

        if m=='chi2':
            values = values_chi2
            errors = errors_chi2
            limits = limits_chi2
        elif m=='bb':
            values = values_bb
            errors = errors_bb
            limits = limits_bb
        elif m=='hb':
            values = values_hb
            errors = errors_hb
            limits = limits_hb
            
        ## Correlations
        ##--------------
        c = []
        
        ## Grids
        c.append([1,3]) # I3.3/I11.3 - I7.7/I11.3
        ## 3.4
        c.append([7,8]) # I3.4/I3.3 - [NeIII]/[NeII]
        # c.append([7,12]) # I3.4/I3.3 - [NeIII]/[NeII](r)

        
        # c.append([2,1]) # I6.2/I11.3 - I3.3/I11.3
        c.append([2,3]) # I6.2/I11.3 - I7.7/I11.3
        c.append([4,3]) # I8.6/I11.3 - I7.7/I11.3

        # c.append([7,15]) # I3.4/I3.3 - [H2S1-7]/7.7+8.6
        # c.append([7,16]) # I3.4/I3.3 - [H2S1]/11.3
        # c.append([1,6]) # I3.3/I11.3 - I17.0/I11.3
        # c.append([1,16]) # I3.3/I11.3 - [H2S1]/11.3
        
        # c.append([7,9]) # I3.4/I3.3 - [SIV]/[SIII]
        # c.append([7,10]) # I3.4/I3.3 - [SIV]/[NeII]
        # c.append([7,13]) # I3.4/I3.3 - [ArIII]/[ArII]

        ## Plot
        ##------
        Ncorr = len(c)
        Nsamp = len(values[0])
        print('Number of correlations: '+str(Ncorr))
        print('Sample size:            '+str(Nsamp))

        for icorr in range(Ncorr):
            filename = path_fig+'corr_'+fnames[c[icorr][0]]+' - '+fnames[c[icorr][1]]+'_'+m+'.png'
                
            ## Major plane
            p = pplot(xlabel=labels[c[icorr][1]], ylabel=labels[c[icorr][0]],
                      xlog=1, ylog=1,
                      # xlim=limits[c[icorr][1]], ylim=limits[c[icorr][0]],
                      # xlim=(2e-3,1e0), ylim=(2e-3,1e0),
                      title=None,
                      figsize=(10,8), left=.12, right=.8, top=.95, bottom=.1,
                      # figsize=(16,8), right=.65, left=.08, bottom=.1, top=.95,
                      titlesize=20, labelsize=20, ticksize=20)

            ## Sample points
            ##---------------
            for i in range(Nsamp):
                cx = values[c[icorr][1]][i]
                cy = values[c[icorr][0]][i]
                if m=='chi2':
                    cxerr = errors[c[icorr][1]][i]
                    cyerr = errors[c[icorr][0]][i]
                else:
                    cxerr = np.array(errors[c[icorr][1]])[:,i][:,np.newaxis]
                    cyerr = np.array(errors[c[icorr][0]])[:,i][:,np.newaxis]

                cy = (cy-.26)*1990
                cyerr = (cyerr+.16)*1990
                p.add_plot(cx, cy, yerr=cyerr, xerr=cxerr,
                           fmt='s', ec='k', zorder=100,
                           marker=mlist[i], markersize=10, capsize=2,
                           c=clist[i], label=spec_name[i,0])

                ## Manage legend
                ph = p.ax.get_legend_handles_labels()[0]
                if i<7:
                    fl = p.ax.legend(handles=ph,
                                     loc='upper left', bbox_to_anchor=(1,1),
                                     fontsize=20, framealpha=0)
                elif i<14:
                    p.ax.add_artist(fl)
                    fl = p.ax.legend(handles=ph[7:],
                                     loc='upper left', bbox_to_anchor=(1.2,1),
                                     fontsize=20, framealpha=0)
                elif i<21:
                    p.ax.add_artist(fl)
                    fl = p.ax.legend(handles=ph[14:],
                                     loc='upper left', bbox_to_anchor=(1.4,1),
                                     fontsize=20, framealpha=0)
                elif i<25:
                    p.ax.add_artist(fl)
                    fl = p.ax.legend(handles=ph[21:],
                                     loc='upper left', bbox_to_anchor=(1,0.54),
                                     fontsize=20, framealpha=0)
                elif i<27:
                    p.ax.add_artist(fl)
                    fl = p.ax.legend(handles=ph[25:],
                                     loc='upper left', bbox_to_anchor=(1.2,0.54),
                                     fontsize=20, framealpha=0)
                elif i<29:
                    p.ax.add_artist(fl)
                    fl = p.ax.legend(handles=ph[27:],
                                     loc='upper left', bbox_to_anchor=(1.4,0.54),
                                     fontsize=20, framealpha=0)
                elif i<31:
                    p.ax.add_artist(fl)
                    fl = p.ax.legend(handles=ph[29:],
                                     loc='upper left', bbox_to_anchor=(1,0.25),
                                     fontsize=20, framealpha=0)
                elif i<33:
                    p.ax.add_artist(fl)
                    fl = p.ax.legend(handles=ph[31:],
                                     loc='upper left', bbox_to_anchor=(1.2,0.25),
                                     fontsize=20, framealpha=0)
                elif i<35:
                    p.ax.add_artist(fl)
                    fl = p.ax.legend(handles=ph[33:],
                                     loc='upper left', bbox_to_anchor=(1.4,0.25),
                                     fontsize=20, framealpha=0)
                elif i<37:
                    p.ax.add_artist(fl)
                    fl = p.ax.legend(handles=ph[35:],
                                     loc='upper left', bbox_to_anchor=(1,0.1),
                                     fontsize=20, framealpha=0)
                elif i<39:
                    p.ax.add_artist(fl)
                    fl = p.ax.legend(handles=ph[37:],
                                     loc='upper left', bbox_to_anchor=(1.2,0.1),
                                     fontsize=20, framealpha=0)
                else:
                    p.ax.add_artist(fl)
                    fl = p.ax.legend(handles=ph[39:],
                                     loc='upper left', bbox_to_anchor=(1.4,0.1),
                                     fontsize=20, framealpha=0)

            # p.ax.xaxis.set_major_formatter(ScalarFormatter()) # major
            p.ax.xaxis.set_minor_formatter(NullFormatter()) # minor
            # p.ax.yaxis.set_major_formatter(ScalarFormatter()) # major
            p.ax.yaxis.set_minor_formatter(NullFormatter()) # minor
            # p.ax.legend(handles=phandles,
            #             loc='upper left', bbox_to_anchor=(1,1),
            #             fontsize=20, framealpha=0)

            ## Grids
            ##-------
            if icorr==0:
                p.ax.set_xlim((1e-1,1e2))
                p.ax.set_ylim((1e-3,1e1))

                ## Ncmin
                grid_Ncmin = [19,99,460]
                color_Ncmin = ['k', 'grey', 'lightgrey']
                for j in range(len(grid_Ncmin)):
                    xgrid = []
                    ygrid = []
                    for igrid in range(Ngrid):
                        if Ncmin[igrid]==grid_Ncmin[j] and ISRFtype[igrid]=='Mathis83':
                        # if Ncmin[igrid]==grid_Ncmin[j] and ISRFtype[igrid]=='BB_T30000K':
                            xgrid.append(I77_I113[igrid])
                            ygrid.append(I33_I113[igrid])
                    xgrid = np.array(xgrid)
                    ygrid = np.array(ygrid)
                    ind = np.argsort(ygrid)
                    p.add_plot(xgrid[ind], ygrid[ind], c=color_Ncmin[j], lw=5,
                               label=r'$N_C^min=$'+str(grid_Ncmin[j]))

                ## Manage legend
                ph = p.ax.get_legend_handles_labels()[0]
                p.ax.add_artist(fl)
                fl = p.ax.legend(handles=ph[41:],
                                 loc='upper left', bbox_to_anchor=(.6,.25),
                                 fontsize=20, framealpha=0)
                
                ## fion
                grid_fion = [0,.3,.5,.7,1]
                # color_fion = ['k', 'dimgrey', 'darkgrey', 'grey', 'silver', 'lightgrey']
                color_fion = ['k', 'grey', 'darkgrey', 'lightgrey', 'whitesmoke']
                for j in range(len(grid_fion)):
                    xgrid = []
                    ygrid = []
                    for igrid in range(Ngrid):
                        if fion[igrid]==grid_fion[j] and ISRFtype[igrid]=='Mathis83':
                        # if fion[igrid]==grid_fion[j] and ISRFtype[igrid]=='BB_T30000K':
                            xgrid.append(I77_I113[igrid])
                            ygrid.append(I33_I113[igrid])
                    xgrid = np.array(xgrid)
                    ygrid = np.array(ygrid)
                    ind = np.argsort(ygrid)
                    p.add_plot(xgrid[ind], ygrid[ind], c=color_fion[j], ls='--', lw=3,
                               label=r'$f^+=$'+str(grid_fion[j]))

                ## Manage legend
                ph = p.ax.get_legend_handles_labels()[0]
                p.ax.add_artist(fl)
                fl = p.ax.legend(handles=ph[44:],
                                 loc='upper left', bbox_to_anchor=(.01,.99),
                                 fontsize=20, framealpha=0)

                ## Arrows
                p.ax.annotate('', xytext=(5e-1,1e-2),xy=(2e0,2e-3),
                              arrowprops=dict(facecolor='k',shrink=0.))
                p.ax.text(3e-1,2e-3,'Charge',fontsize=20)
                
                p.ax.annotate('', xytext=(5e1,1e-1),xy=(2e1,1e-2),
                              arrowprops=dict(facecolor='k',shrink=0.))
                p.ax.text(4e1,2e-2,'Size',fontsize=20)

            elif icorr==1:
                p.ax.set_xlim((1e-2,1e0))
                p.ax.set_ylim((1e-2,1e0))

                p.ax.annotate('', xytext=(0.3,0.1),xy=(0.8,0.1),
                              arrowprops=dict(facecolor='k',shrink=0.),
                              xycoords=p.ax.transAxes)
                p.ax.text(0.4,0.05,'ISRF hardness',fontsize=20,
                          transform=p.ax.transAxes)
                
                p.ax.annotate('', xytext=(0.9,0.7),xy=(0.9,0.2),
                              arrowprops=dict(facecolor='k',shrink=0.),
                              xycoords=p.ax.transAxes)
                p.ax.text(0.93,0.3,'Dehydrogenation', fontsize=20,
                          rotation=90, transform=p.ax.transAxes)

            elif icorr==2:
                p.ax.set_xlim((1 ,1e1))
                p.ax.set_ylim((1e2,1e4))
                p.ax.set_ylabel(r'$G_0/(n_e/1\ cm^{-3})\times \sqrt{T_{gas}/10^3 K}$')
                
                cx = values[c[icorr][1]]
                cy = (values[c[icorr][0]]-.26)*1990
                if m=='chi2':
                    cxerr = errors[c[icorr][1]]
                    cyerr = (errors[c[icorr][0]]+.16)*1990
                else:
                    cxerr = np.mean(np.array(errors[c[icorr][1]]),axis=0)
                    cyerr = (np.mean(np.array(errors[c[icorr][0]]),axis=0)+.16)*1990
                print(cxerr.shape, cx.shape)
                popt, pcov = curve_fit(func, cx, cy,)
                                       # sigma=cyerr)
                perr = np.sqrt(np.diag(pcov))
                xgrid = np.logspace(0,1,10000)
                yy = []
                for i in range(1000):
                    a = np.random.normal(popt[0],perr[0])
                    b = np.random.normal(popt[1],perr[1])
                    yy.append(func(xgrid,a,b))
                yy = np.array(yy)
                fmin = np.nanmin(yy,axis=0)
                fmax = np.nanmax(yy,axis=0)
                
                p.add_plot(xgrid, func(xgrid,*popt), zorder=100,
                           c='white',label='f={:.2f}*x{:.2f}'.format(*popt))
                p.ax.fill_between(xgrid, fmin,fmax, facecolor='grey')
                
                ## Manage legend
                ph = p.ax.get_legend_handles_labels()[0]
                p.ax.add_artist(fl)
                fl = p.ax.legend(handles=ph[41:],
                                 loc='upper left', bbox_to_anchor=(.01,.99),
                                 fontsize=20, framealpha=0)

                
            
            p.save(filename, transparent=True)
            

print('>>> Coucou show_corr [Done] <<<')
