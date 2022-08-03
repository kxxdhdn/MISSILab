#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

## astylo
from astylo.iolib import read_fits, write_fits
from astylo.plib import plot2d_m

## Local
from param_irc import (
    src, fits_irc, path_tests
)
colorlib = ['r', 'orange', 'y', 'g', 'c', 'b', \
	'm', 'pink', 'chocolate', 'lime'] # 10 colors
label = ['1', '2', '3', '4', '5', '6']
path_fig = path_tests+'Y12spec/'
print(fits_irc)
x = read_fits(fits_irc[0]).wave

data0 = read_fits(fits_irc[0]).data
data1 = read_fits(fits_irc[1]).data
data2 = read_fits(fits_irc[2]).data
data3 = read_fits(fits_irc[3]).data
data4 = read_fits(fits_irc[4]).data
data5 = read_fits(fits_irc[5]).data

y0 = [data0[:,0,0], data0[:,12,0]]
plot2d_m(None, y0, x, xlim=[2.5,4.5], ylim=[0, 2.5], 
	cl=colorlib, lablist=label, legend='upper left',
         title='H').save(path_fig+'H')
y1 = [data1[:,0,0], data1[:,4,0], data1[:,8,0],
      data1[:,12,0], data1[:,16,0], data1[:,20,0]]
plot2d_m(None, y1, x, xlim=[2.5,4.5], ylim=[0, 500], 
	cl=colorlib, lablist=label, legend='upper left',
         title='A').save(path_fig+'A')
y2 = [data2[:,0,0], data2[:,12,0]]
plot2d_m(None, y2, x, xlim=[2.5,4.5], ylim=[0, 2.5], 
	cl=colorlib, lablist=label, legend='upper left',
         title='I').save(path_fig+'I')
y3 = [data3[:,0,0], data3[:,4,0], data3[:,8,0],
      data3[:,12,0], data3[:,16,0], data3[:,20,0]]
plot2d_m(None, y3, x, xlim=[2.5,4.5], ylim=[0, 200], 
	cl=colorlib, lablist=label, legend='upper left',
         title='C').save(path_fig+'C')
y4 = [data4[:,0,0], data4[:,12,0]]
plot2d_m(None, y4, x, xlim=[2.5,4.5], ylim=[0, 2.5], 
	cl=colorlib, lablist=label, legend='upper left',
         title='G').save(path_fig+'G')
y5 = [data5[:,0,0], data5[:,4,0], data5[:,8,0],
      data5[:,12,0], data5[:,16,0], data5[:,20,0]]
plot2d_m(None, y5, x, xlim=[2.5,4.5], ylim=[0, 150], 
	cl=colorlib, lablist=label, legend='upper left',
         title='B').save(path_fig+'B')
