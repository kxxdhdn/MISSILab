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

## Settings
mode = ['bb','hb'] # bb befor hb
bb_ok = False
hb_ok = False
locext = 5e-2 # inside
# locext = 5e1 # outside
## lnRband1 (3.3) = 18 (simu)
ihyp = 0 # Python index (+1 for the function "ij2icorr")
## lnRband11 (11.2) = 28 (simu)
ihyp2 = 1 # Python index (+1 for the function "ij2icorr")
## Upper triangle
if ihyp2<ihyp:
    ihyp, ihyp2 = ihyp2, ihyp

Nbin = 'sqrt' # doane
is_rel = input('Reliable/unreliable estimation (y/n)? ')

xlab = 'Lag'
for m in mode:
    if m=='bb':
        MODE = 'non-HB'
        c_acf = 'k'
    elif m=='hb':
        MODE = 'HB'
        c_acf = 'b'
    filout = path_out+'fit_'+m

    h5_out = pathlib.Path(filout+h5ext)
    if h5_out.exists():
        if m=='bb':
            bb_ok = True
        elif m=='hb':
            hb_ok = True

        hyparname = read_hdf5(filout, 'Hyperparameter label')
        Nhypar = len(hyparname)

        ## Parameter names
        parname2 = [hyparname[ihyp],hyparname[ihyp2]]
        for i2, ih in enumerate([ihyp,ihyp2]):
            parname2[i2] = plotname(hyparname[ih])

        ## Label names
        muname0 = r'$\mu$( '+parname2[0]+' )'
        signame0 = r'$\sigma$( '+parname2[0]+' )'
        corrname0 = 'corr( '+parname2[0]+','+parname2[1]+' )'
        ## icorr = (i-1)*N-(i+1)C2+j
        icorr = ij2icorr(ihyp+1, ihyp2+1, Nhypar) # -1 for Python index

## Mu
##----
for m in mode:
    if m=='bb':
        MODE = 'non-HB'
        c_acf = 'k'
    elif m=='hb':
        MODE = 'HB'
        c_acf = 'b'
    filout = path_out+'fit_'+m

    h5_out = pathlib.Path(filout+h5ext)
    if h5_out.exists():

        ## Read ACF analysis
        acf_mu = read_hdf5(filout,
            'Autocorrelation function for mu('+hyparname[ihyp]+')')
        t_int_mu = read_hdf5(filout,
            'Integrated autocorrelation time for the means')
        Neff_mu = read_hdf5(filout,
            'Effective sample size for the means')

        Nlag = len(acf_mu)
        lag = np.arange(Nlag)
        
        ## Plot ACF
        xmin, xmax = 0, Nlag
        ylab = 'ACF of '+muname0
        if m=='bb' or (m=='hb' and not bb_ok):
            p = pplot(lag, acf_mu, c=c_acf, alpha=.5, lw=2, label=MODE,
                      xlog=True, nonposx='sym', xtkform='log_sci',
                      xlabel=xlab, ylabel=ylab, xlim=(xmin,xmax),
                      figsize=(12,4), xysize=20, tksize=20)
            p.ax.axhline(0, c='grey', ls='dashdot', lw=2,)
        else:
            p.add_plot(lag, acf_mu, c=c_acf, alpha=.5, lw=2, label=MODE)
        p.ax.axvline(abs(t_int_mu[ihyp]), c=c_acf, alpha=.5, ls='dashed', lw=2,
                     label=r'$\tau_{\rm int}=$'+f'{abs(t_int_mu[ihyp]):.1f}, ' \
                           +r'$N_{\rm eff}=$'+f'{abs(Neff_mu[ihyp])}')# ('+MODE+')')
        p.append_handles()
        p.set_legend(loc='extra upper right',locext=locext,fontsize=20,framealpha=0)
        
        filename = path_fig+'ACF_'+hyparname[ihyp]+'_mu.png'
        p.save(filename, transparent=True, figtight=True, close=False)
        if (m=='bb' and not hb_ok) or m=='hb':
            p.close()


## Mu (burn-in)
##--------------
for m in mode:
    if m=='bb':
        MODE = 'non-HB'
        c_acf = 'k'
    elif m=='hb':
        MODE = 'HB'
        c_acf = 'b'
    filout = path_out+'fit_'+m

    h5_out = pathlib.Path(filout+h5ext)
    if h5_out.exists():
        
        ## Read ACF analysis
        acf_mu = read_hdf5(filout,
            'Autocorrelation function for mu('+hyparname[ihyp]+') after burn-in')
        t_int_mu = read_hdf5(filout,
            'Integrated autocorrelation time for the means after burn-in')
        Neff_mu = read_hdf5(filout,
            'Effective sample size for the means after burn-in')
        
        Nlag = len(acf_mu)
        lag = np.arange(Nlag)
        
        ## Plot ACF
        xmin, xmax = 0, Nlag
        ylab = 'ACF of '+muname0
        if m=='bb' or (m=='hb' and not bb_ok):
            p = pplot(lag, acf_mu, c=c_acf, alpha=.5, lw=2, label=MODE,
                      xlog=True, nonposx='sym', xtkform='log_sci',
                      xlabel=xlab, ylabel=ylab, xlim=(xmin,xmax),
                      figsize=(12,4), xysize=20, tksize=20)
            p.ax.axhline(0, c='grey', ls='dashdot', lw=2,)
        else:
            p.add_plot(lag, acf_mu, c=c_acf, alpha=.5, lw=2, label=MODE)
        p.ax.axvline(abs(t_int_mu[ihyp]), c=c_acf, alpha=.5, ls='dashed', lw=2,
                      label=r'$\tau_{\rm int}=$'+f'{abs(t_int_mu[ihyp]):.1f}, ' \
                            +r'$N_{\rm eff}=$'+f'{abs(Neff_mu[ihyp])}')# ('+MODE+')')
        p.append_handles()
        p.set_legend(loc='extra upper right',locext=locext,fontsize=20,framealpha=0)
        
        filename = path_fig+'ACF_'+hyparname[ihyp]+'_mu_burn.png'
        p.save(filename, transparent=True, figtight=True, close=False)
        if (m=='bb' and not hb_ok) or m=='hb':
            p.close()


## Sig
##-----
for m in mode:
    if m=='bb':
        MODE = 'non-HB'
        c_acf = 'k'
    elif m=='hb':
        MODE = 'HB'
        c_acf = 'b'
    filout = path_out+'fit_'+m

    h5_out = pathlib.Path(filout+h5ext)
    if h5_out.exists():

        ## Read ACF analysis
        acf_sig = read_hdf5(filout,
            'Autocorrelation function for sig('+hyparname[ihyp]+')')
        t_int_sig = read_hdf5(filout,
            'Integrated autocorrelation time for the standard deviations')
        Neff_sig = read_hdf5(filout,
            'Effective sample size for the standard deviations')

        Nlag = len(acf_sig)
        lag = np.arange(Nlag)
        
        ## Plot ACF
        xmin, xmax = 0, Nlag
        ylab = 'ACF of '+signame0
        if m=='bb' or (m=='hb' and not bb_ok):
            p = pplot(lag, acf_sig, c=c_acf, alpha=.5, lw=2, label=MODE,
                      xlog=True, nonposx='sym', xtkform='log_sci',
                      xlabel=xlab, ylabel=ylab, xlim=(xmin,xmax),
                      figsize=(12,4), xysize=20, tksize=20)
            p.ax.axhline(0, c='grey', ls='dashdot', lw=2,)
        else:
            p.add_plot(lag, acf_sig, c=c_acf, alpha=.5, lw=2, label=MODE)
        p.ax.axvline(abs(t_int_sig[ihyp]), c=c_acf, alpha=.5, ls='dashed', lw=2,
                     label=r'$\tau_{\rm int}=$'+f'{abs(t_int_sig[ihyp]):.1f}, ' \
                           +r'$N_{\rm eff}=$'+f'{abs(Neff_sig[ihyp])}')# ('+MODE+')')
        p.append_handles()
        p.set_legend(loc='extra upper right',locext=locext,fontsize=20,framealpha=0)
        
        filename = path_fig+'ACF_'+hyparname[ihyp]+'_sig.png'
        p.save(filename, transparent=True, figtight=True, close=False)
        if (m=='bb' and not hb_ok) or m=='hb':
            p.close()


## Sig (burn-in)
##---------------
for m in mode:
    if m=='bb':
        MODE = 'non-HB'
        c_acf = 'k'
    elif m=='hb':
        MODE = 'HB'
        c_acf = 'b'
    filout = path_out+'fit_'+m

    h5_out = pathlib.Path(filout+h5ext)
    if h5_out.exists():
        
        ## Read ACF analysis
        acf_sig = read_hdf5(filout,
            'Autocorrelation function for sig('+hyparname[ihyp]+') after burn-in')
        t_int_sig = read_hdf5(filout,
            'Integrated autocorrelation time for the standard deviations after burn-in')
        Neff_sig = read_hdf5(filout,
            'Effective sample size for the standard deviations after burn-in')
        
        Nlag = len(acf_sig)
        lag = np.arange(Nlag)
        
        ## Plot ACF
        xmin, xmax = 0, Nlag
        ylab = 'ACF of '+signame0
        if m=='bb' or (m=='hb' and not bb_ok):
            p = pplot(lag, acf_sig, c=c_acf, alpha=.5, lw=2, label=MODE,
                      xlog=True, nonposx='sym', xtkform='log_sci',
                      xlabel=xlab, ylabel=ylab, xlim=(xmin,xmax),
                      figsize=(12,4), xysize=20, tksize=20)
            p.ax.axhline(0, c='grey', ls='dashdot', lw=2,)
        else:
            p.add_plot(lag, acf_sig, c=c_acf, alpha=.5, lw=2, label=MODE)
        p.ax.axvline(abs(t_int_sig[ihyp]), c=c_acf, alpha=.5, ls='dashed', lw=2,
                      label=r'$\tau_{\rm int}=$'+f'{abs(t_int_sig[ihyp]):.1f}, ' \
                            +r'$N_{\rm eff}=$'+f'{abs(Neff_sig[ihyp])}')# ('+MODE+')')
        p.append_handles()
        p.set_legend(loc='extra upper right',locext=locext,fontsize=20,framealpha=0)
        
        filename = path_fig+'ACF_'+hyparname[ihyp]+'_sig_burn.png'
        p.save(filename, transparent=True, figtight=True, close=False)
        if (m=='bb' and not hb_ok) or m=='hb':
            p.close()


## Corr
##------
for m in mode:
    if m=='bb':
        MODE = 'non-HB'
        c_acf = 'k'
    elif m=='hb':
        MODE = 'HB'
        c_acf = 'b'
    filout = path_out+'fit_'+m

    h5_out = pathlib.Path(filout+h5ext)
    if h5_out.exists() and icorr>0:

        ## Read ACF analysis
        acf_corr = read_hdf5(filout,
            'Autocorrelation function for corr('+str(icorr)+')')
        t_int_corr = read_hdf5(filout,
            'Integrated autocorrelation time for the correlations')
        Neff_corr = read_hdf5(filout,
            'Effective sample size for the correlations')

        Nlag = len(acf_corr)
        lag = np.arange(Nlag)
        
        ## Plot ACF
        xmin, xmax = 0, Nlag
        ylab = 'ACF of '+corrname0
        if m=='bb' or (m=='hb' and not bb_ok):
            p = pplot(lag, acf_corr, c=c_acf, alpha=.5, lw=2, label=MODE,
                      xlog=True, nonposx='sym', xtkform='log_sci',
                      xlabel=xlab, ylabel=ylab, xlim=(xmin,xmax),
                      figsize=(12,4), xysize=20, tksize=20)
            p.ax.axhline(0, c='grey', ls='dashdot', lw=2,)
        else:
            p.add_plot(lag, acf_corr, c=c_acf, alpha=.5, lw=2, label=MODE)
        p.ax.axvline(abs(t_int_corr[icorr]), c=c_acf, alpha=.5, ls='dashed', lw=2,
                     label=r'$\tau_{\rm int}=$'+f'{t_int_corr[icorr]:.1f}, ' \
                           +r'$N_{\rm eff}=$'+f'{Neff_corr[icorr]}')# ('+MODE+')')
        p.append_handles()
        p.set_legend(loc='extra upper right',locext=locext,fontsize=20,framealpha=0)
        
        filename = path_fig+'ACF_'+hyparname[ihyp]+'_corr'+str(icorr)+'.png'
        p.save(filename, transparent=True, figtight=True, close=False)
        if (m=='bb' and not hb_ok) or m=='hb':
            p.close()


## Corr (burn-in)
##----------------
for m in mode:
    if m=='bb':
        MODE = 'non-HB'
        c_acf = 'k'
    elif m=='hb':
        MODE = 'HB'
        c_acf = 'b'
    filout = path_out+'fit_'+m

    h5_out = pathlib.Path(filout+h5ext)
    if h5_out.exists() and icorr>0:
        
        ## Read ACF analysis
        acf_corr = read_hdf5(filout,
            'Autocorrelation function for corr('+str(icorr)+') after burn-in')
        t_int_corr = read_hdf5(filout,
            'Integrated autocorrelation time for the correlations after burn-in')
        Neff_corr = read_hdf5(filout,
            'Effective sample size for the correlations after burn-in')
        
        Nlag = len(acf_corr)
        lag = np.arange(Nlag)
        
        ## Plot ACF
        xmin, xmax = 0, Nlag
        ylab = 'ACF of '+corrname0
        if m=='bb' or (m=='hb' and not bb_ok):
            p = pplot(lag, acf_corr, c=c_acf, alpha=.5, lw=2, label=MODE,
                      xlog=True, nonposx='sym', xtkform='log_sci',
                      xlabel=xlab, ylabel=ylab, xlim=(xmin,xmax),
                      figsize=(12,4), xysize=20, tksize=20)
            p.ax.axhline(0, c='grey', ls='dashdot', lw=2,)
        else:
            p.add_plot(lag, acf_corr, c=c_acf, alpha=.5, lw=2, label=MODE)
        p.ax.axvline(abs(t_int_corr[icorr]), c=c_acf, alpha=.5, ls='dashed', lw=2,
                      label=r'$\tau_{\rm int}=$'+f'{abs(t_int_corr[icorr]):.1f}, ' \
                            +r'$N_{\rm eff}=$'+f'{abs(Neff_corr[icorr])}')# ('+MODE+')')
        p.append_handles()
        p.set_legend(loc='extra upper right',locext=locext,fontsize=20,framealpha=0)

        filename = path_fig+'ACF_'+hyparname[ihyp]+'_corr'+str(icorr)+'_burn.png'
        p.save(filename, transparent=True, figtight=True, close=False)
        if (m=='bb' and not hb_ok) or m=='hb':
            p.close()


## tau_int distribution
##----------------------
for m in mode:
    if m=='bb':
        MODE = 'non-HB'
        c_acf = 'k'
    elif m=='hb':
        MODE = 'HB'
        c_acf = 'b'
    filout = path_out+'fit_'+m

    h5_out = pathlib.Path(filout+h5ext)
    if h5_out.exists() and (is_rel=='y' or is_rel=='n'):

        ## Mu
        t_int_mu = read_hdf5(filout,
            'Integrated autocorrelation time for the means after burn-in')
        if is_rel=='y':
            mask = t_int_mu>0 # reliable
            rel = 'reliable'
        else:
            mask = np.full(len(t_int_mu), True, dtype=bool) # unreliable
            rel = 'unreliable'
        xmin = min(abs(t_int_mu[mask]))
        xmax = max(abs(t_int_mu[mask]))
        p = plotool(1, 1, figsize=(10,8))
        p.set_ax([1,0],xlabel=r'$\tau_{\rm int}$ of $\mu$',ylabel='Counts',
                 xlim=(xmin,xmax),# xlog=1, xtkform='mylog', ytkform='log_sci',
                 xsize=20,ysize=20,xtksize=20,ytksize=20,)
        count, bins, patches = p.ax.hist(abs(t_int_mu[mask]), bins=Nbin, log=1,
            facecolor='c', edgecolor='None', alpha=.5, label=r'$\mu$ (After burn-in)')
        p.append_handles()
        p.set_legend(loc='extra upper right',locext=1e-1,fontsize=20,framealpha=0)

        filename = path_fig+'acf_tint_mu_'+m+'_'+rel+'.png'
        p.save(filename, transparent=True, figtight=True)
        
        ## Sig
        t_int_sig = read_hdf5(filout,
            'Integrated autocorrelation time for the standard deviations after burn-in')
        if is_rel=='y':
            mask = t_int_sig>0 # reliable
            rel = 'reliable'
        else:
            mask = np.full(len(t_int_sig), True, dtype=bool) # unreliable
            rel = 'unreliable'
        xmin = min(abs(t_int_sig[mask]))
        xmax = max(abs(t_int_sig[mask]))
        p = plotool(1, 1, figsize=(10,8))
        p.set_ax([1,0],xlabel=r'$\tau_{\rm int}$ of $\sigma$',ylabel='Counts',
                 xlim=(xmin,xmax),# xlog=1, xtkform='mylog', ytkform='log_sci',
                 xsize=20,ysize=20,xtksize=20,ytksize=20,)
        count, bins, patches = p.ax.hist(abs(t_int_sig[mask]), bins=Nbin, log=1,
            facecolor='c', edgecolor='None', alpha=.5, label=r'$\sigma$ (After burn-in)')
        p.append_handles()
        p.set_legend(loc='extra upper right',locext=1e-1,fontsize=20,framealpha=0)

        filename = path_fig+'acf_tint_sig_'+m+'_'+rel+'.png'
        p.save(filename, transparent=True, figtight=True)
        
        ## Corr
        t_int_corr = read_hdf5(filout,
            'Integrated autocorrelation time for the correlations after burn-in')
        if is_rel=='y':
            mask = t_int_corr>0 # reliable
            rel = 'reliable'
        else:
            mask = np.full(len(t_int_corr), True, dtype=bool) # unreliable
            rel = 'unreliable'
        xmin = min(abs(t_int_corr[mask]))
        xmax = max(abs(t_int_corr[mask]))
        p = plotool(1, 1, figsize=(10,8))
        p.set_ax([1,0],xlabel=r'$\tau_{\rm int}$ of $\rho$',ylabel='Counts',
                 xlim=(xmin,xmax),# xlog=1, xtkform='mylog', ytkform='log_sci',
                 xsize=20,ysize=20,xtksize=20,ytksize=20,)
        count, bins, patches = p.ax.hist(abs(t_int_corr[mask]), bins=Nbin, log=1,
            facecolor='c', edgecolor='None', alpha=.5,
            label=r'Correlation coefficients, $\rho$ (After burn-in)')
        p.append_handles()
        p.set_legend(loc='extra upper right',locext=1e-1,fontsize=20,framealpha=0)

        filename = path_fig+'acf_tint_corr_'+m+'_'+rel+'.png'
        p.save(filename, transparent=True, figtight=True)


print('>>> Coucou show_acf [Done] <<<')
