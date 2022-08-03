#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This is the visualization of correlations

"""

import os, pathlib
import numpy as np
import matplotlib.colors as mcolors

## rapyuta
from rapyuta.arrays import closest
from rapyuta.inout import read_hdf5

## local
from auxil import croot


## Path
##------
path_out = croot+'../out/'
h5_obs = path_out+'observation_MIR'
h5_model = path_out+'input_model'

## Read fit
##----------
labB = read_hdf5(h5_model, 'label band')
labL = read_hdf5(h5_model, 'label line')

## Colors and markers
##--------------------
spec_name = read_hdf5(h5_obs, 'spectrum labels')
clist = []
mlist = []
for y in range(spec_name.shape[0]):
    if spec_name[y,0][0]=='A':
        clist.append( ((y)/14+.5, 0, 0) )
        mlist.append('*')
    elif spec_name[y,0][0]=='B':
        clist.append( (.8, (y-7+1)/14+.2, .2) )
        mlist.append('*')
    elif spec_name[y,0][0]=='C':
        clist.append( ((y-14+1)/14+.5, (y-14+1)/14+.5, 0) )
        mlist.append('*')
    elif spec_name[y,0][0]=='D':
        clist.append( (0, (y-21)/8+.5, 0) )
        mlist.append('d')
    # elif spec_name[y,0][0]=='E':
    #     clist.append('lime')
    #     mlist.append('d')
    elif spec_name[y,0][0]=='F':
        clist.append( (.9, (y-25+1)/6+.5, .8) )
        mlist.append('v')
    elif spec_name[y,0][0]=='G':
        clist.append( ((y-27+1)/3+.3, 0, (y-27+1)/3+.3) )
        mlist.append('v')
    elif spec_name[y,0][0]=='H':
        clist.append( (0, (y-29+1)/3+.3, (y-29+1)/3+.3) )
        mlist.append('v')
    elif spec_name[y,0][0]=='I':
        clist.append( (0, 0, (y-31+1)/3+.3) )
        mlist.append('v')
    elif spec_name[y,0][0]=='J':
        clist.append( (0, 0, (y-33+1)/3+.3) )
        mlist.append('^')
    elif spec_name[y,0][0]=='K':
        clist.append( ((y-35+1)/3+.3, 0, (y-35+1)/3+.3) )
        mlist.append('^')
    elif spec_name[y,0][0]=='L':
        clist.append( (.9, (y-37+1)/6+.5, .8) )
        mlist.append('^')
    elif spec_name[y,0][0]=='M':
        clist.append( (0, (y-39+1)/3+.3, (y-39+1)/3+.3) )
        mlist.append('^')
    # elif spec_name[y,0][0]=='N':
    #     clist.append('b')
    #     mlist.append('s')

## Ratios
##--------
mode = ['chi2', 'bb', 'hb']

for m in mode:
    filout = path_out+'fit_'+m
    h5_out = pathlib.Path(filout+'.h5')
    if h5_out.exists():

        labels = []
        fnames = []
        
        if m=='chi2':
            parname = read_hdf5(filout, 'Parameter label')
            par = read_hdf5(filout, 'Best fitted parameter value')
            parerr = read_hdf5(filout, 'Best fitted parameter error')
            par[par==0] = np.nan
            
            values_chi2 = []
            errors_chi2 = []
            limits_chi2 = []
        else:
            filmcmc = path_out+'parlog_fit_'+m
            parname = read_hdf5(filmcmc, 'Parameter label')
            parmcmc = read_hdf5(filmcmc, 'Parameter values')
            t_end = read_hdf5(filmcmc, 'Last index')[0]
            t_burnin = int(t_end/10) - 1
            
            if m=='bb':
                values_bb = []
                errors_bb = []
                limits_bb = []
            else:
                values_hb = []
                errors_hb = []
                limits_hb = []

        ## 11.3 complex
        i1 = np.where(labB=='Main 11.2')[0][0]+1
        ipar1 = np.where(parname=='lnRband'+str(i1))[0][0]
        i2 = np.where(labB=='Plateau 11.3')[0][0]+1
        ipar2 = np.where(parname=='lnRband'+str(i2))[0][0]
        if m=='chi2':
            v0 = 1 + np.exp(par[ipar2,:,:])
            e0 = (np.exp(par[ipar2,:,:]) * parerr[ipar2,:,:])**2
        else:
            v0 = 1 + np.exp(parmcmc[t_burnin:t_end,ipar2,:,:])

        ##----------------------------------------------

        ## I11.3/I11.2
        if m=='chi2':
            v = v0.flatten()
            e = np.sqrt( e0 ).flatten()
            values_chi2.append(v)
            errors_chi2.append(e)
            limits_chi2.append((1e-2,1e1))
        else:
            vf = v0
            v = np.median( vf, axis=0 ).flatten()
            e1 = v - np.quantile( vf, 1/3, axis=0 ).flatten()
            e2 = np.quantile( vf, 2/3, axis=0 ).flatten() - v
            e = [e1,e2]
            if m=='bb':
                values_bb.append(v)
                errors_bb.append(e)
                limits_bb.append((1e-2,1e1))
            else:
                values_hb.append(v)
                errors_hb.append(e)
                limits_hb.append((1e-2,1e1))
        labels.append('I11.3/I11.2')
        fnames.append('I11.3_ov_I11.2')


        ## I3.3/I11.3
        i1 = np.where(labB=='Main 3.3')[0][0]+1
        ipar1 = np.where(parname=='lnRband'+str(i1))[0][0]
        if m=='chi2':
            v = ( np.exp(par[ipar1,:,:]) / v0 ).flatten()
            e = ( np.sqrt( (np.exp(par[ipar1,:,:])*parerr[ipar1,:,:])**2
                           + np.exp(par[ipar1,:,:])**2 * e0 )
                  / v0 ).flatten()
            values_chi2.append(v)
            errors_chi2.append(e)
            limits_chi2.append((1e-2,1e1))
        else:
            vf = np.exp(parmcmc[t_burnin:t_end,ipar1,:,:]) / v0
            v = np.median( vf, axis=0 ).flatten()
            e1 = v - np.quantile( vf, 1/3, axis=0 ).flatten()
            e2 = np.quantile( vf, 2/3, axis=0 ).flatten() - v
            e = [e1,e2]
            if m=='bb':
                values_bb.append(v)
                errors_bb.append(e)
                limits_bb.append((1e-2,1e1))
            else:
                values_hb.append(v)
                errors_hb.append(e)
                limits_hb.append((1e-2,1e1))
        labels.append('I3.3/I11.3')
        fnames.append('I3.3_ov_I11.3')


        ## I6.2/I11.3
        i1 = np.where(labB=='Main 6.2 (1)')[0][0]+1
        ipar1 = np.where(parname=='lnRband'+str(i1))[0][0]
        i2 = np.where(labB=='Main 6.2 (2)')[0][0]+1
        ipar2 = np.where(parname=='lnRband'+str(i2))[0][0]
        if m=='chi2':
            v = ( (np.exp(par[ipar1,:,:])+np.exp(par[ipar2,:,:])) / v0 ).flatten()
            e = ( np.sqrt( (np.exp(par[ipar1,:,:])*parerr[ipar1,:,:])**2
                           + (np.exp(par[ipar2,:,:])*parerr[ipar2,:,:])**2
                           + (np.exp(par[ipar1,:,:])+np.exp(par[ipar2,:,:]))**2 * e0 )
                  / v0 ).flatten()
            values_chi2.append(v)
            errors_chi2.append(e)
            limits_chi2.append((1e-2,1e1))
        else:
            vf = np.exp(parmcmc[t_burnin:t_end,ipar1,:,:]
                        +parmcmc[t_burnin:t_end,ipar2,:,:]) / v0
            v = np.median( vf, axis=0 ).flatten()
            e1 = v - np.quantile( vf, 1/3, axis=0 ).flatten()
            e2 = np.quantile( vf, 2/3, axis=0 ).flatten() - v
            e = [e1,e2]
            if m=='bb':
                values_bb.append(v)
                errors_bb.append(e)
                limits_bb.append((1e-2,1e1))
            else:
                values_hb.append(v)
                errors_hb.append(e)
                limits_hb.append((1e-2,1e1))
        labels.append('I6.2/I11.3')
        fnames.append('I6.2_ov_I11.3')
            
        ## I7.7/I11.3
        i1 = np.where(labB=='Plateau 7.7')[0][0]+1
        ipar1 = np.where(parname=='lnRband'+str(i1))[0][0]
        i2 = np.where(labB=='Main 7.7 (1)')[0][0]+1
        ipar2 = np.where(parname=='lnRband'+str(i2))[0][0]
        i3 = np.where(labB=='Main 7.7 (2)')[0][0]+1
        ipar3 = np.where(parname=='lnRband'+str(i3))[0][0]
        if m=='chi2':
            v = ( (np.exp(par[ipar1,:,:])+np.exp(par[ipar2,:,:])+np.exp(par[ipar3,:,:])) / v0 ).flatten()
            e = ( np.sqrt( (np.exp(par[ipar1,:,:])*parerr[ipar1,:,:])**2
                           + (np.exp(par[ipar2,:,:])*parerr[ipar2,:,:])**2
                           + (np.exp(par[ipar3,:,:])*parerr[ipar3,:,:])**2
                           + (np.exp(par[ipar1,:,:])+np.exp(par[ipar2,:,:])+np.exp(par[ipar3,:,:]))**2 * e0 )
                  / v0 ).flatten()
            values_chi2.append(v)
            errors_chi2.append(e)
            limits_chi2.append((1e-2,1e1))
        else:
            vf = np.exp(parmcmc[t_burnin:t_end,ipar1,:,:]
                        +parmcmc[t_burnin:t_end,ipar2,:,:]
                        +parmcmc[t_burnin:t_end,ipar3,:,:]) / v0
            v = np.median( vf, axis=0 ).flatten()
            e1 = v - np.quantile( vf, 1/3, axis=0 ).flatten()
            e2 = np.quantile( vf, 2/3, axis=0 ).flatten() - v
            e = [e1,e2]
            if m=='bb':
                values_bb.append(v)
                errors_bb.append(e)
                limits_bb.append((1e-2,1e1))
            else:
                values_hb.append(v)
                errors_hb.append(e)
                limits_hb.append((1e-2,1e1))
        labels.append('I7.7/I11.3')
        fnames.append('I7.7_ov_I11.3')

        ## I8.6/I11.3
        i1 = np.where(labB=='Main 8.6')[0][0]+1
        ipar1 = np.where(parname=='lnRband'+str(i1))[0][0]
        if m=='chi2':
            v = ( np.exp(par[ipar1,:,:]) / v0 ).flatten()
            e = ( np.sqrt( (np.exp(par[ipar1,:,:])*parerr[ipar1,:,:])**2
                           + np.exp(par[ipar1,:,:])**2 * e0 )
                  / v0 ).flatten()
            values_chi2.append(v)
            errors_chi2.append(e)
            limits_chi2.append((1e-2,1e1))
        else:
            vf = np.exp(parmcmc[t_burnin:t_end,ipar1,:,:]) / v0
            v = np.median( vf, axis=0 ).flatten()
            e1 = v - np.quantile( vf, 1/3, axis=0 ).flatten()
            e2 = np.quantile( vf, 2/3, axis=0 ).flatten() - v
            e = [e1,e2]
            if m=='bb':
                values_bb.append(v)
                errors_bb.append(e)
                limits_bb.append((1e-2,1e1))
            else:
                values_hb.append(v)
                errors_hb.append(e)
                limits_hb.append((1e-2,1e1))
        labels.append('I8.6/I11.3')
        fnames.append('I8.6_ov_I11.3')

        ## I12.7/I11.3
        i1 = np.where(labB=='Main 12.7 (1)')[0][0]+1
        ipar1 = np.where(parname=='lnRband'+str(i1))[0][0]
        i2 = np.where(labB=='Main 12.7 (2)')[0][0]+1
        ipar2 = np.where(parname=='lnRband'+str(i2))[0][0]
        if m=='chi2':
            v = ( (np.exp(par[ipar1,:,:])+np.exp(par[ipar2,:,:])) / v0 ).flatten()
            e = ( np.sqrt( (np.exp(par[ipar1,:,:])*parerr[ipar1,:,:])**2
                           + (np.exp(par[ipar2,:,:])*parerr[ipar2,:,:])**2
                           + (np.exp(par[ipar1,:,:])+np.exp(par[ipar2,:,:]))**2 * e0 )
                  / v0 ).flatten()
            values_chi2.append(v)
            errors_chi2.append(e)
            limits_chi2.append((1e-2,1e1))
        else:
            vf = np.exp(parmcmc[t_burnin:t_end,ipar1,:,:]
                        +parmcmc[t_burnin:t_end,ipar2,:,:]) / v0
            v = np.median( vf, axis=0 ).flatten()
            e1 = v - np.quantile( vf, 1/3, axis=0 ).flatten()
            e2 = np.quantile( vf, 2/3, axis=0 ).flatten() - v
            e = [e1,e2]
            if m=='bb':
                values_bb.append(v)
                errors_bb.append(e)
                limits_bb.append((1e-2,1e1))
            else:
                values_hb.append(v)
                errors_hb.append(e)
                limits_hb.append((1e-2,1e1))
        labels.append('I12.7/I11.3')
        fnames.append('I12.7_ov_I11.3')

        ## I17.0/I11.3
        i1 = np.where(labB=='Plateau 17.0')[0][0]+1
        ipar1 = np.where(parname=='lnRband'+str(i1))[0][0]
        if m=='chi2':
            v = ( np.exp(par[ipar1,:,:]) / v0 ).flatten()
            e = ( np.sqrt( (np.exp(par[ipar1,:,:])*parerr[ipar1,:,:])**2
                           + np.exp(par[ipar1,:,:])**2 * e0 )
                  / v0 ).flatten()
            values_chi2.append(v)
            errors_chi2.append(e)
            limits_chi2.append((1e-2,1e1))
        else:
            vf = np.exp(parmcmc[t_burnin:t_end,ipar1,:,:]) / v0
            v = np.median( vf, axis=0 ).flatten()
            e1 = v - np.quantile( vf, 1/3, axis=0 ).flatten()
            e2 = np.quantile( vf, 2/3, axis=0 ).flatten() - v
            e = [e1,e2]
            if m=='bb':
                values_bb.append(v)
                errors_bb.append(e)
                limits_bb.append((1e-2,1e1))
            else:
                values_hb.append(v)
                errors_hb.append(e)
                limits_hb.append((1e-2,1e1))
        labels.append('I17.0/I11.3')
        fnames.append('I17.0_ov_I11.3')
            
        ## I3.4/I3.3
        i1 = np.where(labB=='Main 3.4')[0][0]+1
        ipar1 = np.where(parname=='lnRband'+str(i1))[0][0]
        i2 = np.where(labB=='Main 3.3')[0][0]+1
        ipar2 = np.where(parname=='lnRband'+str(i2))[0][0]
        if m=='chi2':
            v = np.exp( par[ipar1,:,:]-par[ipar2,:,:] ).flatten()
            e = v * np.sqrt( parerr[ipar1,:,:]**2+parerr[ipar2,:,:]**2 ).flatten()
            values_chi2.append(v)
            errors_chi2.append(e)
            limits_chi2.append((1e-2,1e1))
        else:
            vf = np.exp(parmcmc[t_burnin:t_end,ipar1,:,:]
                        -parmcmc[t_burnin:t_end,ipar2,:,:])
            v = np.median( vf, axis=0 ).flatten()
            e1 = v - np.quantile( vf, 1/3, axis=0 ).flatten()
            e2 = np.quantile( vf, 2/3, axis=0 ).flatten() - v
            e = [e1,e2]
            if m=='bb':
                values_bb.append(v)
                errors_bb.append(e)
                limits_bb.append((1e-2,1e1))
            else:
                values_hb.append(v)
                errors_hb.append(e)
                limits_hb.append((1e-2,1e1))
        labels.append('I3.4/I3.3')
        fnames.append('I3.4_ov_I3.3')

        ## [NeIII]/[NeII]
        i1 = np.where(labL=='NeIII1')[0][0]+1
        ipar1 = np.where(parname=='lnRline'+str(i1))[0][0]
        i2 = np.where(labL=='NeII  ')[0][0]+1
        ipar2 = np.where(parname=='lnRline'+str(i2))[0][0]
        if m=='chi2':
            v = np.exp( par[ipar1,:,:]-par[ipar2,:,:] ).flatten()
            e = v * np.sqrt( parerr[ipar1,:,:]**2+parerr[ipar2,:,:]**2 ).flatten()
            values_chi2.append(v)
            errors_chi2.append(e)
            limits_chi2.append((1e-2,1e1))
        else:
            vf = np.exp(parmcmc[t_burnin:t_end,ipar1,:,:]
                        -parmcmc[t_burnin:t_end,ipar2,:,:])
            v = np.median( vf, axis=0 ).flatten()
            e1 = v - np.quantile( vf, 1/3, axis=0 ).flatten()
            e2 = np.quantile( vf, 2/3, axis=0 ).flatten() - v
            e = [e1,e2]
            if m=='bb':
                values_bb.append(v)
                errors_bb.append(e)
                limits_bb.append((1e-2,1e1))
            else:
                values_hb.append(v)
                errors_hb.append(e)
                limits_hb.append((1e-2,1e1))
        labels.append('[NeIII]/[NeII]')
        fnames.append('[NeIII]_ov_[NeII]')

        ## [SIV]/[SIII]
        i1 = np.where(labL=='SIV   ')[0][0]+1
        ipar1 = np.where(parname=='lnRline'+str(i1))[0][0]
        i2 = np.where(labL=='SIII1 ')[0][0]+1
        ipar2 = np.where(parname=='lnRline'+str(i2))[0][0]
        if m=='chi2':
            v = np.exp( par[ipar1,:,:]-par[ipar2,:,:] ).flatten()
            e = v * np.sqrt( parerr[ipar1,:,:]**2+parerr[ipar2,:,:]**2 ).flatten()
            values_chi2.append(v)
            errors_chi2.append(e)
            limits_chi2.append((1e-2,1e1))
        else:
            vf = np.exp(parmcmc[t_burnin:t_end,ipar1,:,:]
                        -parmcmc[t_burnin:t_end,ipar2,:,:])
            v = np.median( vf, axis=0 ).flatten()
            e1 = v - np.quantile( vf, 1/3, axis=0 ).flatten()
            e2 = np.quantile( vf, 2/3, axis=0 ).flatten() - v
            e = [e1,e2]
            if m=='bb':
                values_bb.append(v)
                errors_bb.append(e)
                limits_bb.append((1e-2,1e1))
            else:
                values_hb.append(v)
                errors_hb.append(e)
                limits_hb.append((1e-2,1e1))
        labels.append('[SIV]/[SIII]')
        fnames.append('[SIV]_ov_[SIII]')

        ## [SIV]/[NeII]
        i1 = np.where(labL=='SIV   ')[0][0]+1
        ipar1 = np.where(parname=='lnRline'+str(i1))[0][0]
        i2 = np.where(labL=='NeII  ')[0][0]+1
        ipar2 = np.where(parname=='lnRline'+str(i2))[0][0]
        if m=='chi2':
            v = np.exp( par[ipar1,:,:]-par[ipar2,:,:] ).flatten()
            e = v * np.sqrt( parerr[ipar1,:,:]**2+parerr[ipar2,:,:]**2 ).flatten()
            values_chi2.append(v)
            errors_chi2.append(e)
            limits_chi2.append((1e-2,1e1))
        else:
            vf = np.exp(parmcmc[t_burnin:t_end,ipar1,:,:]
                        -parmcmc[t_burnin:t_end,ipar2,:,:])
            v = np.median( vf, axis=0 ).flatten()
            e1 = v - np.quantile( vf, 1/3, axis=0 ).flatten()
            e2 = np.quantile( vf, 2/3, axis=0 ).flatten() - v
            e = [e1,e2]
            if m=='bb':
                values_bb.append(v)
                errors_bb.append(e)
                limits_bb.append((1e-2,1e1))
            else:
                values_hb.append(v)
                errors_hb.append(e)
                limits_hb.append((1e-2,1e1))
        labels.append('[SIV]/[NeII]')
        fnames.append('[SIV]_ov_[NeII]')
            
        ## [NeIII]/[SIII]
        i1 = np.where(labL=='NeIII1')[0][0]+1
        ipar1 = np.where(parname=='lnRline'+str(i1))[0][0]
        i2 = np.where(labL=='SIII1 ')[0][0]+1
        ipar2 = np.where(parname=='lnRline'+str(i2))[0][0]
        if m=='chi2':
            v = np.exp( par[ipar1,:,:]-par[ipar2,:,:] ).flatten()
            e = v * np.sqrt( parerr[ipar1,:,:]**2+parerr[ipar2,:,:]**2 ).flatten()
            values_chi2.append(v)
            errors_chi2.append(e)
            limits_chi2.append((1e-2,1e1))
        else:
            vf = np.exp(parmcmc[t_burnin:t_end,ipar1,:,:]
                        -parmcmc[t_burnin:t_end,ipar2,:,:])
            v = np.median( vf, axis=0 ).flatten()
            e1 = v - np.quantile( vf, 1/3, axis=0 ).flatten()
            e2 = np.quantile( vf, 2/3, axis=0 ).flatten() - v
            e = [e1,e2]
            if m=='bb':
                values_bb.append(v)
                errors_bb.append(e)
                limits_bb.append((1e-2,1e1))
            else:
                values_hb.append(v)
                errors_hb.append(e)
                limits_hb.append((1e-2,1e1))
        labels.append('[NeIII]/[SIII]')
        fnames.append('[NeIII]_ov_[SIII]')

        ## [NeIII]/[NeII](r)
        i1 = np.where(labL=='NeIII1')[0][0]+1
        ipar1 = np.where(parname=='lnRline'+str(i1))[0][0]
        i2 = np.where(labL=='NeII  ')[0][0]+1
        ipar2 = np.where(parname=='lnRline'+str(i2))[0][0]
        i3 = np.where(labB=='Main 12.7 (1)')[0][0]+1
        ipar3 = np.where(parname=='lnRband'+str(i3))[0][0]
        if m=='chi2':
            v = np.exp( par[ipar1,:,:]-par[ipar2,:,:]+par[ipar3,:,:] ).flatten()
            e = v * np.sqrt( parerr[ipar1,:,:]**2+parerr[ipar2,:,:]**2+parerr[ipar3,:,:]**2 ).flatten()
            values_chi2.append(v)
            errors_chi2.append(e)
            limits_chi2.append((1e-2,1e1))
        else:
            vf = np.exp(parmcmc[t_burnin:t_end,ipar1,:,:]
                        -parmcmc[t_burnin:t_end,ipar2,:,:]
                        +parmcmc[t_burnin:t_end,ipar3,:,:])
            v = np.median( vf, axis=0 ).flatten()
            e1 = v - np.quantile( vf, 1/3, axis=0 ).flatten()
            e2 = np.quantile( vf, 2/3, axis=0 ).flatten() - v
            e = [e1,e2]
            if m=='bb':
                values_bb.append(v)
                errors_bb.append(e)
                limits_bb.append((1e-2,1e1))
            else:
                values_hb.append(v)
                errors_hb.append(e)
                limits_hb.append((1e-2,1e1))
        labels.append('[NeIII]/[NeII](r)')
        fnames.append('[NeIII]_ov_[NeII](r)')

        ## [ArIII]/[ArII]
        i1 = np.where(labL=='ArIII1')[0][0]+1
        ipar1 = np.where(parname=='lnRline'+str(i1))[0][0]
        i2 = np.where(labL=='ArII  ')[0][0]+1
        ipar2 = np.where(parname=='lnRline'+str(i2))[0][0]
        if m=='chi2':
            v = np.exp( par[ipar1,:,:]-par[ipar2,:,:] ).flatten()
            e = v * np.sqrt( parerr[ipar1,:,:]**2+parerr[ipar2,:,:]**2 ).flatten()
            values_chi2.append(v)
            errors_chi2.append(e)
            limits_chi2.append((1e-2,1e1))
        else:
            vf = np.exp(parmcmc[t_burnin:t_end,ipar1,:,:]
                        -parmcmc[t_burnin:t_end,ipar2,:,:])
            v = np.median( vf, axis=0 ).flatten()
            e1 = v - np.quantile( vf, 1/3, axis=0 ).flatten()
            e2 = np.quantile( vf, 2/3, axis=0 ).flatten() - v
            e = [e1,e2]
            if m=='bb':
                values_bb.append(v)
                errors_bb.append(e)
                limits_bb.append((1e-2,1e1))
            else:
                values_hb.append(v)
                errors_hb.append(e)
                limits_hb.append((1e-2,1e1))
        labels.append('[ArIII]/[ArII]')
        fnames.append('[ArIII]_ov_[ArII]')

        ## [SIV]/[NeII](r)
        i1 = np.where(labL=='SIV   ')[0][0]+1
        ipar1 = np.where(parname=='lnRline'+str(i1))[0][0]
        i2 = np.where(labL=='NeII  ')[0][0]+1
        ipar2 = np.where(parname=='lnRline'+str(i2))[0][0]
        i3 = np.where(labB=='Main 12.7 (1)')[0][0]+1
        ipar3 = np.where(parname=='lnRband'+str(i3))[0][0]
        if m=='chi2':
            v = np.exp( par[ipar1,:,:]-par[ipar2,:,:]+par[ipar3,:,:] ).flatten()
            e = v * np.sqrt( parerr[ipar1,:,:]**2+parerr[ipar2,:,:]**2+parerr[ipar3,:,:]**2 ).flatten()
            values_chi2.append(v)
            errors_chi2.append(e)
            limits_chi2.append((1e-2,1e1))
        else:
            vf = np.exp(parmcmc[t_burnin:t_end,ipar1,:,:]
                        -parmcmc[t_burnin:t_end,ipar2,:,:]
                        +parmcmc[t_burnin:t_end,ipar3,:,:])
            v = np.median( vf, axis=0 ).flatten()
            e1 = v - np.quantile( vf, 1/3, axis=0 ).flatten()
            e2 = np.quantile( vf, 2/3, axis=0 ).flatten() - v
            e = [e1,e2]
            if m=='bb':
                values_bb.append(v)
                errors_bb.append(e)
                limits_bb.append((1e-2,1e1))
            else:
                values_hb.append(v)
                errors_hb.append(e)
                limits_hb.append((1e-2,1e1))
        labels.append('[SIV]/[NeII](r)')
        fnames.append('[SIV]_ov_[NeII](r)')

        ## [H2S1-7]/(I7.7+I8.6)
        i1 = np.where(labL=='H2S1  ')[0][0]+1
        ipar1 = np.where(parname=='lnRline'+str(i1))[0][0]
        i2 = np.where(labL=='H2S2  ')[0][0]+1
        ipar2 = np.where(parname=='lnRline'+str(i2))[0][0]
        i3 = np.where(labL=='H2S3  ')[0][0]+1
        ipar3 = np.where(parname=='lnRline'+str(i3))[0][0]
        i4 = np.where(labL=='H2S5  ')[0][0]+1
        ipar4 = np.where(parname=='lnRline'+str(i4))[0][0]
        i5 = np.where(labL=='H2S7  ')[0][0]+1
        ipar5 = np.where(parname=='lnRline'+str(i5))[0][0]
        i6 = np.where(labB=='Plateau 7.7')[0][0]+1
        ipar6 = np.where(parname=='lnRband'+str(i6))[0][0]
        i7 = np.where(labB=='Main 7.7 (1)')[0][0]+1
        ipar7 = np.where(parname=='lnRband'+str(i7))[0][0]
        i8 = np.where(labB=='Main 7.7 (2)')[0][0]+1
        ipar8 = np.where(parname=='lnRband'+str(i8))[0][0]
        i9 = np.where(labB=='Main 8.6')[0][0]+1
        ipar9 = np.where(parname=='lnRband'+str(i9))[0][0]
        if m=='chi2':
            v = np.exp( par[ipar1,:,:]+par[ipar2,:,:]+par[ipar3,:,:]
                       +par[ipar4,:,:]+par[ipar5,:,:]-par[ipar6,:,:]
                       -par[ipar7,:,:]-par[ipar8,:,:]-par[ipar9,:,:]).flatten()
            e = v * np.sqrt( parerr[ipar1,:,:]**2+parerr[ipar2,:,:]**2+parerr[ipar3,:,:]**2
                            +parerr[ipar4,:,:]**2+parerr[ipar5,:,:]**2+parerr[ipar6,:,:]**2
                            +parerr[ipar7,:,:]**2+parerr[ipar8,:,:]**2+parerr[ipar9,:,:]**2).flatten()
            values_chi2.append(v)
            errors_chi2.append(e)
            limits_chi2.append((1e-2,1e1))
        else:
            vf = np.exp( parmcmc[t_burnin:t_end,ipar1,:,:]
                        +parmcmc[t_burnin:t_end,ipar2,:,:]
                        +parmcmc[t_burnin:t_end,ipar3,:,:]
                        +parmcmc[t_burnin:t_end,ipar4,:,:]
                        +parmcmc[t_burnin:t_end,ipar5,:,:]
                        -parmcmc[t_burnin:t_end,ipar6,:,:]
                        -parmcmc[t_burnin:t_end,ipar7,:,:]
                        -parmcmc[t_burnin:t_end,ipar8,:,:]
                        -parmcmc[t_burnin:t_end,ipar9,:,:])
            v = np.median( vf, axis=0 ).flatten()
            e1 = v - np.quantile( vf, 1/3, axis=0 ).flatten()
            e2 = np.quantile( vf, 2/3, axis=0 ).flatten() - v
            e = [e1,e2]
            if m=='bb':
                values_bb.append(v)
                errors_bb.append(e)
                limits_bb.append((1e-2,1e1))
            else:
                values_hb.append(v)
                errors_hb.append(e)
                limits_hb.append((1e-2,1e1))
        labels.append('[H2S1-7]/(I7.7+I8.6)')
        fnames.append('[H2S1-7]_ov_I7.7+I8.6')

        ## [H2S1]/I11.3
        i1 = np.where(labL=='H2S1  ')[0][0]+1
        ipar1 = np.where(parname=='lnRline'+str(i1))[0][0]
        if m=='chi2':
            v = ( np.exp(par[ipar1,:,:]) / v0 ).flatten()
            e = ( np.sqrt( (np.exp(par[ipar1,:,:])*parerr[ipar1,:,:])**2
                           + np.exp(par[ipar1,:,:])**2 * e0 )
                  / v0 ).flatten()
            values_chi2.append(v)
            errors_chi2.append(e)
            limits_chi2.append((1e-2,1e1))
        else:
            vf = np.exp(parmcmc[t_burnin:t_end,ipar1,:,:]) / v0
            v = np.median( vf, axis=0 ).flatten()
            e1 = v - np.quantile( vf, 1/3, axis=0 ).flatten()
            e2 = np.quantile( vf, 2/3, axis=0 ).flatten() - v
            e = [e1,e2]
            if m=='bb':
                values_bb.append(v)
                errors_bb.append(e)
                limits_bb.append((1e-2,1e1))
            else:
                values_hb.append(v)
                errors_hb.append(e)
                limits_hb.append((1e-2,1e1))
        labels.append('[H2S1]/I11.3')
        fnames.append('[H2S1]_ov_I11.3')
        

print(len(values_chi2))
print(len(values_hb))
