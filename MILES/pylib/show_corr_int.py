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


snames = ['A1-7',
          'B1-7',
          'C1-7',
          'D1-4',
          'F1-2',
          'G1-2',
          'H1-2',
          'I1-2',
          'J1-2',
          'K1-2',
          'L1-2',
          'M1-2',]

clist = []
mlist = []
for y in range(len(snames)):
    if y==0:
        clist.append('r') # A
        mlist.append('*')
    elif y==1:
        clist.append('orange')
        mlist.append('*')
    elif y==2:
        clist.append('y')
        mlist.append('*')
    elif y==3:
        clist.append('g') # D
        mlist.append('d')
    elif y==4:
        clist.append('pink') # F
        mlist.append('v')
    elif y==5:
        clist.append('m')
        mlist.append('v')
    elif y==6:
        clist.append('c')
        mlist.append('v')
    elif y==7:
        clist.append('b') # I
        mlist.append('v')
    elif y==8:
        clist.append('b')
        mlist.append('^')
    elif y==9:
        clist.append('m') # K
        mlist.append('^')
    elif y==10:
        clist.append('pink')
        mlist.append('^')
    elif y==11:
        clist.append('c')
        mlist.append('^')




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
        
        c.append([7,8]) # I3.4/I3.3 - [NeIII]/[NeII]
        c.append([7,15]) # I3.4/I3.3 - [H2S1-7]
        # c.append([7,12]) # I3.4/I3.3 - [NeIII]/[NeII](r)
        # c.append([7,9]) # I3.4/I3.3 - [SIV]/[SIII]
        # c.append([7,10]) # I3.4/I3.3 - [SIV]/[NeII]
        # c.append([7,13]) # I3.4/I3.3 - [ArIII]/[ArII]
        # c.append([7,14]) # I3.4/I3.3 - [ArIII]/[ArII]

        ## Plot
        ##------
        Ncorr = len(c)
        Nsamp = len(values[0])
        print('Number of correlations: '+str(Ncorr))
        print('Sample size:            '+str(Nsamp))

        for icorr in range(Ncorr):
            filename = path_fig+'cint_corr_'+fnames[c[icorr][0]]+' - '+fnames[c[icorr][1]]+'_'+m+'.png'
                
            ## Major plane
            p = pplot(xlabel=labels[c[icorr][1]], ylabel=labels[c[icorr][0]],
                      xlog=1, ylog=1,
                      # xlim=limits[c[icorr][1]], ylim=limits[c[icorr][0]],
                      # xlim=(1e-2,1e0), ylim=(1e-2,1e0),
                      # xlim=(1e-4,1e2), ylim=(1e-4,1e2),
                      title=None,
                      figsize=(10,8), left=.12, right=.8, top=.95, bottom=.1,
                      titlesize=20, labelsize=20, ticksize=20)

            ## Sample points
            ##---------------
            Ax, Ay, Axerr, Ayerr = 0, 0, 0, 0
            Bx, By, Bxerr, Byerr = 0, 0, 0, 0
            Cx, Cy, Cxerr, Cyerr = 0, 0, 0, 0
            Dx, Dy, Dxerr, Dyerr = 0, 0, 0, 0
            Fx, Fy, Fxerr, Fyerr = 0, 0, 0, 0
            Gx, Gy, Gxerr, Gyerr = 0, 0, 0, 0
            Hx, Hy, Hxerr, Hyerr = 0, 0, 0, 0
            Ix, Iy, Ixerr, Iyerr = 0, 0, 0, 0
            Jx, Jy, Jxerr, Jyerr = 0, 0, 0, 0
            Kx, Ky, Kxerr, Kyerr = 0, 0, 0, 0
            Lx, Ly, Lxerr, Lyerr = 0, 0, 0, 0
            Mx, My, Mxerr, Myerr = 0, 0, 0, 0
            for i in range(Nsamp):
                cx = values[c[icorr][1]][i]
                cy = values[c[icorr][0]][i]
                if m=='chi2':
                    cxerr = errors[c[icorr][1]][i]
                    cyerr = errors[c[icorr][0]][i]
                else:
                    cxerr = np.array(errors[c[icorr][1]])[:,i][:,np.newaxis]
                    cyerr = np.array(errors[c[icorr][0]])[:,i][:,np.newaxis]

                ux = 1e40
                dx = 1e-40
                uy = 1e40
                dy = 1e-40
                if i<7 and cx>dx and cx<ux and cy>dy and cy<uy:
                    if i==0:
                        countA = 1
                    else:
                        countA += 1
                    Ax += cx
                    Axerr += cxerr
                    Ay += cy
                    Ayerr += cyerr
                elif i<14 and cx>dx and cx<ux and cy>dy and cy<uy:
                    if i==7:
                        countB = 1
                    else:
                        countB += 1
                    Bx += cx
                    Bxerr += cxerr
                    By += cy
                    Byerr += cyerr
                elif i<21 and cx>dx and cx<ux and cy>dy and cy<uy:
                    if i==14:
                        countC = 1
                    else:
                        countC += 1
                    Cx += cx
                    Cxerr += cxerr
                    Cy += cy
                    Cyerr += cyerr
                elif i<25 and cx>dx and cx<ux and cy>dy and cy<uy:
                    if i==21:
                        countD = 1
                    else:
                        countD += 1
                    Dx += cx
                    Dxerr += cxerr
                    Dy += cy
                    Dyerr += cyerr
                elif i<27 and cx>dx and cx<ux and cy>dy and cy<uy:
                    if i==25:
                        countF = 1
                    else:
                        countF += 1
                    Fx += cx
                    Fxerr += cxerr
                    Fy += cy
                    Fyerr += cyerr
                elif i<29 and cx>dx and cx<ux and cy>dy and cy<uy:
                    if i==27:
                        countG = 1
                    else:
                        countG += 1
                    Gx += cx
                    Gxerr += cxerr
                    Gy += cy
                    Gyerr += cyerr
                elif i<31 and cx>dx and cx<ux and cy>dy and cy<uy:
                    if i==29:
                        countH = 1
                    else:
                        countH += 1
                    Hx += cx
                    Hxerr += cxerr
                    Hy += cy
                    Hyerr += cyerr
                elif i<33 and cx>dx and cx<ux and cy>dy and cy<uy:
                    if i==31:
                        countI = 1
                    else:
                        countI += 1
                    Ix += cx
                    Ixerr += cxerr
                    Iy += cy
                    Iyerr += cyerr
                elif i<35 and cx>dx and cx<ux and cy>dy and cy<uy:
                    if i==33:
                        countJ = 1
                    else:
                        countJ += 1
                    Jx += cx
                    Jxerr += cxerr
                    Jy += cy
                    Jyerr += cyerr
                elif i<37 and cx>dx and cx<ux and cy>dy and cy<uy:
                    if i==35:
                        countK = 1
                    else:
                        countK += 1
                    Kx += cx
                    Kxerr += cxerr
                    Ky += cy
                    Kyerr += cyerr
                elif i<39 and cx>dx and cx<ux and cy>dy and cy<uy:
                    if i==37:
                        countL = 1
                    else:
                        countL += 1
                    Lx += cx
                    Lxerr += cxerr
                    Ly += cy
                    Lyerr += cyerr
                elif cx>dx and cx<ux and cy>dy and cy<uy:
                    if i==39:
                        countM = 1
                    else:
                        countM += 1
                    Mx += cx
                    Mxerr += cxerr
                    My += cy
                    Myerr += cyerr

            for j in range (len(snames)):
                if j==0:
                    x    = Ax    /countA
                    xerr = Axerr /countA
                    y    = Ay    /countA
                    yerr = Ayerr /countA
                elif j==1:      
                    x    = Bx    /countB
                    xerr = Bxerr /countB
                    y    = By    /countB
                    yerr = Byerr /countB
                elif j==2:      
                    x    = Cx    /countC
                    xerr = Cxerr /countC
                    y    = Cy    /countC
                    yerr = Cyerr /countC
                elif j==3:      
                    x    = Dx    /countD
                    xerr = Dxerr /countD
                    y    = Dy    /countD
                    yerr = Dyerr /countD
                elif j==4:      
                    x    = Fx    /countF
                    xerr = Fxerr /countF
                    y    = Fy    /countF
                    yerr = Fyerr /countF
                elif j==5:       
                    x    = Gx    /countG
                    xerr = Gxerr /countG
                    y    = Gy    /countG
                    yerr = Gyerr /countG
                elif j==6:       
                    x    = Hx    /countH
                    xerr = Hxerr /countH
                    y    = Hy    /countH
                    yerr = Hyerr /countH
                elif j==7:       
                    x    = Ix    /countI
                    xerr = Ixerr /countI
                    y    = Iy    /countI
                    yerr = Iyerr /countI
                elif j==8:       
                    x    = Jx    /countJ
                    xerr = Jxerr /countJ
                    y    = Jy    /countJ
                    yerr = Jyerr /countJ
                elif j==9:       
                    x    = Kx    /countK
                    xerr = Kxerr /countK
                    y    = Ky    /countK
                    yerr = Kyerr /countK
                elif j==10:      
                    x    = Lx    /countL
                    xerr = Lxerr /countL
                    y    = Ly    /countL
                    yerr = Lyerr /countL
                elif j==11:     
                    x    = Mx    /countM
                    xerr = Mxerr /countM
                    y    = My    /countM
                    yerr = Myerr /countM

                    
                p.add_plot(x, y, yerr=yerr, xerr=xerr,
                           fmt='s', ec='grey', zorder=100,
                           marker=mlist[j], markersize=10, capsize=2,
                           c=clist[j], label=snames[j])
                
            # p.ax.xaxis.set_major_formatter(ScalarFormatter()) # major
            p.ax.xaxis.set_minor_formatter(NullFormatter()) # minor
            # p.ax.yaxis.set_major_formatter(ScalarFormatter()) # major
            p.ax.yaxis.set_minor_formatter(NullFormatter()) # minor
            p.ax.legend(#handles=phandles,
                        loc='upper left', bbox_to_anchor=(1,1),
                        fontsize=20, framealpha=0)

            # p.ax.set_xlim((1e-2,1e0))
            # p.ax.set_ylim((1e-2,1e0))

            # p.ax.annotate('', xytext=(0.3,0.1),xy=(0.8,0.1),
            #               arrowprops=dict(facecolor='k',shrink=0.),
            #               xycoords=p.ax.transAxes)
            # p.ax.text(0.4,0.05,'ISRF hardness',fontsize=20,
            #           transform=p.ax.transAxes)
            
            # p.ax.annotate('', xytext=(0.9,0.7),xy=(0.9,0.2),
            #               arrowprops=dict(facecolor='k',shrink=0.),
            #               xycoords=p.ax.transAxes)
            # p.ax.text(0.93,0.3,'Dehydrogenation',
            #           rotation=90,fontsize=20,
            #           transform=p.ax.transAxes)
                
            
            p.save(filename, transparent=True)
            

print('>>> Coucou show_corr [Done] <<<')
