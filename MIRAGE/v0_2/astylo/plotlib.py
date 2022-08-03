#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Plot library

    plotool:
        figure, set_border, Cartesian2d, set_ax,
        plot, set_font, fill, text, save, show
    plot2D, plot2D_m

"""
from astropy import units as u
import numpy as np
from scipy import optimize
import matplotlib as mpl
import matplotlib.pyplot as plt

global colib, LSList
# cmap = mpl.cm.viridis
# norm = mpl.colors.Normalize(vmin=0, vmax=1)
colib = ['k', 'grey', 'r', 'orange', 'y', 'g', 'c', 'b', \
    'm', 'pink', 'chocolate', 'lime'] # 12 colors
# colib = ['blue']*100
LSList = ['-', '--', '-.', ':']
sizeXS, sizeS, sizeM = 4, 6, 8
sizeL, sizeXL, sizeXXL = 10, 12, 14

class plotool:
    """
    PLOT tOOL (2D Cartesian coordinates)
    """
    def __init__(self, x=np.zeros(2), y=np.zeros(2), \
        xerr=None, yerr=None):
        
        # INPUTS
        self.x = x
        self.y = y
        self.xerr = xerr
        self.yerr = yerr

        self.ax = None

    def figure(self, figsize=None, figint=True, \
        nrows=1, ncols=1, sharex=False, sharey=False):
        
        if figint==True:
            plt.ion()

        self.nrows = nrows
        self.ncols = ncols

        if nrows==1 and ncols==1:
            self.fig, self.ax = plt.subplots(nrows, ncols, \
                sharex, sharey, figsize=figsize)
        else:
            self.fig, self.axes = plt.subplots(nrows, ncols, \
                sharex, sharey, figsize=figsize)
        
    def set_border(self, left=None, bottom=None, \
        right=None, top=None, wspace=None, hspace=None):

        plt.subplots_adjust(left=left, bottom=bottom, \
            right=right, top=top, wspace=wspace, hspace=hspace)

    def Cartesian2d(self, c=None, ls=None, lw=None, \
        ec=None, elw=None, lab=None, legend=None, title=None):

        self.markers, self.caps, self.bars = self.ax.errorbar(
            self.x, self.y, yerr=self.yerr, xerr=self.xerr, c=c, \
            ls=ls, lw=lw, ecolor=ec, elinewidth=elw, label=lab)
        
        self.ax.set_title(title)
        if legend is not None:
            self.ax.legend(loc=legend)

    def set_ax(self, xlog=False, ylog=False, \
        xlim=(None,None), ylim=(None,None), xlab=None, ylab=None):

        if xlog==True:
            self.ax.set_xscale('symlog', nonposx='clip')
        if ylog==True:
            self.ax.set_yscale('symlog', nonposy='clip')
        
        self.ax.set_xlim(xlim[0], xlim[1])
        self.ax.set_ylim(ylim[0], ylim[1])

        # self.ax.set_xticks()
        # self.ax.set_yticks()
        # self.ax.set_xticklabels()
        # self.ax.set_yticklabels()
        
        self.ax.set_xlabel(xlab)
        self.ax.set_ylabel(ylab)
    
    def plot(self, nrow=1, ncol=1, \
        xlim=(None, None), ylim=(None, None), \
        xlog=False, ylog=False, c=None, ls=None, lw=None, \
        ec=None, elw=None, xlab=None, ylab=None, \
        lab=None, legend=None, title=None, mod='car2'):

        if mod=='car2':
            if self.nrows==1 and self.ncols==1:
                self.Cartesian2d(c, ls, lw, ec, elw, lab, legend, title)
            else:
                self.ax = self.axes[nrow-1,ncol-1]
                self.Cartesian2d(c, ls, lw, ec, elw, lab, legend, title)
            
            self.set_ax(xlog, ylog, xlim, ylim, xlab, ylab)
        
        else:
            print('*******************')
            print('...Prochainement...')
            print('*******************')

    def set_font(self, fontsize=sizeM, subtitlesize=sizeM, \
        axesize=sizeS, xticksize=sizeS, yticksize=sizeS, \
        legendsize=sizeM, figtitlesize=sizeL):

        plt.rc('font', size=fontsize)            # general text
        plt.rc('axes', titlesize=subtitlesize)    # axes title
        plt.rc('axes', labelsize=axesize)        # x and y labels
        plt.rc('xtick', labelsize=xticksize)    # x tick
        plt.rc('ytick', labelsize=yticksize)    # y tick
        plt.rc('legend', fontsize=legendsize)    # legend
        plt.rc('figure', titlesize=figtitlesize)# figure title

    def fill(self, edgec, fc, falpha=.5):

        self.ax.fill(self.x, self.y, ec=edgec, fc=fc, alpha=falpha)
        
    def text(self, textin, textx, texty):

        self.ax.text(textx, texty, color='b', s=textin)

    def save(self, savename, transparent=False):

        if savename is not None:
            self.fig.savefig(savename, transparent=transparent)

    def show(self):

        plt.ioff()
        plt.show()

##-----------------------------------------------

##                    sub-classes

##-----------------------------------------------

class plot2D(plotool):
    """
    Single 2D curve plot
    """
    def __init__(self, x, y, xerr=None, yerr=None, \
        xlim=(None, None), ylim=(None,None), xlog=None, ylog=None, \
        c='b', ls=None, lw=.5, ec='r', elw=.8, \
        xlab='X', ylab='Y', lab=None, legend=None, \
        title='2D Curve', figsize=None, figint=False, \
        left=.1, bottom=.1, right=.99, top=.9, \
        wspace=None, hspace=None, \
        dofill=False, edgec='g', fc='g', falpha=.5):
        super().__init__(x, y, xerr, yerr)

        self.figure(figsize, figint)

        self.set_border(left=left, bottom=bottom, \
            right=right, top=top, wspace=wspace, hspace=hspace)

        self.Cartesian2d(c, ls, lw, ec, elw, lab, legend, title)

        self.set_ax(xlog, ylog, xlim, ylim, xlab, ylab)

        self.set_font()

        if dofill==True:
            self.fill(edgec, fc, falpha)

class plot2D_m(plotool):
    """
    Multiple 2D curves plot
    """
    def __init__(self, xlist, ylist, xall=None, \
        xerr=None, yerr=None, xlim=(None, None), \
        ylim=(None, None), xlog=False, ylog=False, \
        cl='default', lslist=None, lablist=None, \
        c='b', ls=None, lw=.5, ec='r', elw=.8, \
        xlab='X', ylab='Y', lab=None, legend=None, \
        title=None, figsize=None, figint=False):
        super().__init__()

        self.figure(figsize, figint)

        for i, self.y in enumerate(ylist):
            if xall is not None:
                self.x = xall
            else:
                self.x = xlist[i]

            if xerr is not None:
                self.xerr = xerr[i]
            if yerr is not None:
                self.yerr = yerr[i]
            
            if cl=='default':
                c = colib[i]
            elif cl is not None:
                c = cl[i]

            if lslist is not None:
                ls = lslist[i]

            if lablist is not None:
                lab = lablist[i]

            self.Cartesian2d(c, ls, lw, ec, elw, lab, legend, title)
            
            self.set_ax(xlog, ylog, xlim, ylim, xlab, ylab)
