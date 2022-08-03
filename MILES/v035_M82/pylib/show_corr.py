#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This is the visualization of correlations

"""

import os, pathlib
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

## laputan
from laputan.inout import read_hdf5
from laputan.plots import pplot

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


mode = ['chi2', 'bb', 'hb']

## Path
##------
path_out = croot+'../out/'
path_fig = path_out+'Figures/'
if not os.path.exists(path_fig):
    os.makedirs(path_fig)
h5_model = path_out+'input_model'
h5_analysis = path_out+'input_analysis'

## Read fit
##----------
labB = read_hdf5(h5_model, 'label band')
labL = read_hdf5(h5_model, 'label line')

for m in mode:
    filout = path_out+'fit_'+m
    h5_out = pathlib.Path(filout+'.h5')
    if h5_out.exists():

        ## Correlations
        ##--------------
        corr_x = []
        corr_y = []
        err_x = []
        err_y = []
        label_x = []
        label_y = []
        corrname = []
        
        ##======
        ## Chi2 (suppose all par non-correlated)
        ##======
        if m=='chi2':
            parname = read_hdf5(filout, 'Parameter label')
            par = read_hdf5(filout, 'Best fitted parameter value')
            parerr = read_hdf5(filout, 'Best fitted parameter error')
            par[par==0] = np.nan
            
            ## I11.3/I3.3 - I7.7/I11.3
            iband = np.where(labB=='Main 3.3')[0][0]+1
            ipar = np.where(parname=='lnRband'+str(iband))[0][0]
            corr_y.append(np.exp(-par[ipar,:,:]).flatten())
            err_y.append(
                non_corr_df1(-par[ipar,:,:],parerr[ipar,:,:]).flatten())
            label_y.append('I11.3/I3.3')
            iband = np.where(labB=='Main 7.7 (1)')[0][0]+1
            ipar = np.where(parname=='lnRband'+str(iband))[0][0]
            corr_x.append(np.exp(par[ipar,:,:]).flatten())
            err_x.append(
                non_corr_df1(par[ipar,:,:],parerr[ipar,:,:]).flatten())
            label_x.append('I7.7/I11.3')
            corrname.append('I11.3/I3.3 - I7.7/I11.3')
            ## I3.4/I3.3 - [SIV]/[NeII]
            iband1 = np.where(labB=='Main 3.4')[0][0]+1
            ipar1 = np.where(parname=='lnRband'+str(iband1))[0][0]
            iband2 = np.where(labB=='Main 3.3')[0][0]+1
            ipar2 = np.where(parname=='lnRband'+str(iband2))[0][0]
            corr_y.append(np.exp(
                par[ipar1,:,:]-par[ipar2,:,:]).flatten())
            err_y.append(
                non_corr_df2(par[ipar1,:,:],par[ipar2,:,:],
                             parerr[ipar1,:,:],parerr[ipar2,:,:]).flatten())
            label_y.append('I3.4/I3.3')
            iline1 = np.where(labL=='SIV   ')[0][0]+1
            ipar1 = np.where(parname=='lnRline'+str(iline1))[0][0]
            iline2 = np.where(labL=='NeII  ')[0][0]+1
            ipar2 = np.where(parname=='lnRline'+str(iline2))[0][0]
            corr_x.append(np.exp(
                par[ipar1,:,:]-par[ipar2,:,:]).flatten())
            err_x.append(
                non_corr_df2(par[ipar1,:,:],par[ipar2,:,:],
                             parerr[ipar1,:,:],parerr[ipar2,:,:]).flatten())
            label_x.append('[SIV]/[NeII]')
            corrname.append('I3.4/I3.3 - [SIV]/[NeII]')
            ## I8.6/I11.3 - I6.2/I11.3
            iband = np.where(labB=='Main 8.6')[0][0]+1
            ipar = np.where(parname=='lnRband'+str(iband))[0][0]
            corr_y.append(np.exp(par[ipar,:,:]).flatten())
            err_y.append(
                non_corr_df1(par[ipar,:,:],parerr[ipar,:,:]).flatten())
            label_y.append('I8.6/I11.3')
            iband = np.where(labB=='Main 6.2 (1)')[0][0]+1
            ipar = np.where(parname=='lnRband'+str(iband))[0][0]
            corr_x.append(np.exp(par[ipar,:,:]).flatten())
            err_x.append(
                non_corr_df1(par[ipar,:,:],parerr[ipar,:,:]).flatten())
            label_x.append('I6.2/I11.3')
            corrname.append('I8.6/I11.3 - I6.2/I11.3')
            ## I7.7/I11.3 - I6.2/I11.3
            iband = np.where(labB=='Main 7.7 (1)')[0][0]+1
            ipar = np.where(parname=='lnRband'+str(iband))[0][0]
            corr_y.append(np.exp(par[ipar,:,:]).flatten())
            err_y.append(
                non_corr_df1(par[ipar,:,:],parerr[ipar,:,:]).flatten())
            label_y.append('I7.7/I11.3')
            iband = np.where(labB=='Main 6.2 (1)')[0][0]+1
            ipar = np.where(parname=='lnRband'+str(iband))[0][0]
            corr_x.append(np.exp(par[ipar,:,:]).flatten())
            err_x.append(
                non_corr_df1(par[ipar,:,:],parerr[ipar,:,:]).flatten())
            label_x.append('I6.2/I11.3')
            corrname.append('I7.7/I11.3 - I6.2/I11.3')
            
        ##=======
        ## BB/HB
        ##=======
        else:
            filmcmc = path_out+'parlog_fit_'+m
            parname = read_hdf5(filmcmc, 'Parameter label')
            parmcmc = read_hdf5(filmcmc, 'Parameter values')
            t_end = read_hdf5(filmcmc, 'Last index')[0]
            t_burnin = int(t_end/10) - 1
            
            ## I11.3/I3.3 - I7.7/I11.3
            iband = np.where(labB=='Main 3.3')[0][0]+1
            ipar = np.where(parname=='lnRband'+str(iband))[0][0]
            corr_y.append(np.mean(np.exp(
                -parmcmc[t_burnin:t_end,ipar,:,:]),axis=0).flatten())
            err_y.append(np.std(np.exp(
                -parmcmc[t_burnin:t_end,ipar,:,:]),axis=0).flatten())
            label_y.append('I11.3/I3.3')
            iband = np.where(labB=='Main 7.7 (1)')[0][0]+1
            ipar = np.where(parname=='lnRband'+str(iband))[0][0]
            corr_x.append(np.mean(np.exp(
                parmcmc[t_burnin:t_end,ipar,:,:]),axis=0).flatten())
            err_x.append(np.std(np.exp(
                parmcmc[t_burnin:t_end,ipar,:,:]),axis=0).flatten())
            label_x.append('I7.7/I11.3')
            corrname.append('I11.3/I3.3 - I7.7/I11.3')
            ## I3.4/I3.3 - [SIV]/[NeII]
            iband1 = np.where(labB=='Main 3.4')[0][0]+1
            ipar1 = np.where(parname=='lnRband'+str(iband1))[0][0]
            iband2 = np.where(labB=='Main 3.3')[0][0]+1
            ipar2 = np.where(parname=='lnRband'+str(iband2))[0][0]
            corr_y.append(np.mean(np.exp(
                parmcmc[t_burnin:t_end,ipar1,:,:]-parmcmc[t_burnin:t_end,ipar2,:,:]),
                axis=0).flatten())
            err_y.append(np.std(np.exp(
                parmcmc[t_burnin:t_end,ipar1,:,:]-parmcmc[t_burnin:t_end,ipar2,:,:]),
                axis=0).flatten())
            label_y.append('I3.4/I3.3')
            iline1 = np.where(labL=='SIV   ')[0][0]+1
            ipar1 = np.where(parname=='lnRline'+str(iline1))[0][0]
            iline2 = np.where(labL=='NeII  ')[0][0]+1
            ipar2 = np.where(parname=='lnRline'+str(iline2))[0][0]
            corr_x.append(np.mean(np.exp(
                parmcmc[t_burnin:t_end,ipar1,:,:]-parmcmc[t_burnin:t_end,ipar2,:,:]),
                axis=0).flatten())
            err_x.append(np.std(np.exp(
                parmcmc[t_burnin:t_end,ipar1,:,:]-parmcmc[t_burnin:t_end,ipar2,:,:]),
                axis=0).flatten())
            label_x.append('[SIV]/[NeII]')
            corrname.append('I3.4/I3.3 - [SIV]/[NeII]')
            ## I8.6/I11.3 - I6.2/I11.3
            iband = np.where(labB=='Main 8.6')[0][0]+1
            ipar = np.where(parname=='lnRband'+str(iband))[0][0]
            corr_y.append(np.mean(np.exp(
                parmcmc[t_burnin:t_end,ipar,:,:]),axis=0).flatten())
            err_y.append(np.std(np.exp(
                parmcmc[t_burnin:t_end,ipar,:,:]),axis=0).flatten())
            label_y.append('I8.6/I11.3')
            iband = np.where(labB=='Main 6.2 (1)')[0][0]+1
            ipar = np.where(parname=='lnRband'+str(iband))[0][0]
            corr_x.append(np.mean(np.exp(
                parmcmc[t_burnin:t_end,ipar,:,:]),axis=0).flatten())
            err_x.append(np.std(np.exp(
                parmcmc[t_burnin:t_end,ipar,:,:]),axis=0).flatten())
            label_x.append('I6.2/I11.3')
            corrname.append('I8.6/I11.3 - I6.2/I11.3')
            ## I7.7/I11.3 - I6.2/I11.3
            iband = np.where(labB=='Main 7.7 (1)')[0][0]+1
            ipar = np.where(parname=='lnRband'+str(iband))[0][0]
            corr_y.append(np.mean(np.exp(
                parmcmc[t_burnin:t_end,ipar,:,:]),axis=0).flatten())
            err_y.append(np.std(np.exp(
                parmcmc[t_burnin:t_end,ipar,:,:]),axis=0).flatten())
            label_y.append('I7.7/I11.3')
            iband = np.where(labB=='Main 6.2 (1)')[0][0]+1
            ipar = np.where(parname=='lnRband'+str(iband))[0][0]
            corr_x.append(np.mean(np.exp(
                parmcmc[t_burnin:t_end,ipar,:,:]),axis=0).flatten())
            err_x.append(np.std(np.exp(
                parmcmc[t_burnin:t_end,ipar,:,:]),axis=0).flatten())
            label_x.append('I6.2/I11.3')
            corrname.append('I7.7/I11.3 - I6.2/I11.3')

            # qpar = read_hdf5(filout, 'Quantiles of parameter value')
            # parval = qpar[1]
            # parerrp = qpar[0]
            # parerrn = qpar[2]
            
        Ncorr = np.size(corrname)
        Nsamp = np.size(corr_x[0])

        for icorr in range(Ncorr):
            ## Remove NaNs
            cx = []
            cy = []
            cxe = []
            cye = []
            for ix, x in enumerate(corr_x[icorr]):
                if not (np.isnan(x) or np.isnan(corr_y[icorr][ix])):
                    cx.append(x)
                    cy.append(corr_y[icorr][ix])
                    cxe.append(err_x[icorr][ix])
                    cye.append(err_y[icorr][ix])
            
            filename = path_fig+'corr_'+str(icorr+1)+'_'+m+'.png'

            p = pplot(cx, cy, yerr=cye, xerr=cxe,
                      fmt='o', markersize=.8, c='k', ec='r',
                      xlog=1, ylog=1, nonposx='clip', nonposy='clip',
                      # xlim=(1.e-1,1.e1), ylim=(1.e-1,1.1e2),
                      xlab=label_x[icorr], ylab=label_y[icorr],
                      figsize=(8,8),
                      legend='upper left', label=m, title='M82')

            ## Log linear fit
            popt, pcov = curve_fit(func,cx,cy)#,sigma=cye)
            # print(popt)
            xarr = np.array(sorted(cx))
            
            p.add_plot(xarr, func(xarr,*popt),
                       fmt='y-',
                       label='y={:.2}x+{:.2}'.format(popt[0],popt[1]))
            
            p.save(filename)



# plt.show()

print('>>> Coucou show_corr [Done] <<<')
