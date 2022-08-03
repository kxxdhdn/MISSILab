#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This is the interface of MILES generating following input files:

  /out/set_input.h5       -- input path
  /out/observation_MIR.h5 -- data + initial param
  /out/input_master.h5    -- general param
  /out/input_model.h5     -- model param
  /out/input_extra.h5     -- extra param
  /out/input_analysis.h5  -- analysis param

"""

import os, pathlib
import math
import numpy as np
import subprocess as SP

## rapyuta
from rapyuta.inout import read_fits, write_hdf5, read_hdf5, h5ext

## local
from auxil import (croot, mroot,
                   res, TABLine, TABand, partuning)
         
## Path
##------
h5_path = mroot+'out/set_input'
dirin = croot+'../out/' # define input dir
if not os.path.exists(dirin):
    os.makedirs(dirin)
write_hdf5(h5_path, 'input dir', [dirin], verbose=True)

h5_analysis = dirin+'input_analysis'

program = 'sim_analysis'
noisy = False # verbose/debug for this routine
if noisy==False:
    devnull = open(os.devnull, 'w')
else:
    devnull = None
            
##-----------------------------
##
## Write fit inputs
##
##-----------------------------
input_fit = croot+'input_simu_hb.py'
SP.call('python '+input_fit,
        shell=True, cwd=dirin, stdout=devnull, stderr=SP.STDOUT)


##--------------------------------
##
## Write input_analysis.h5
##
##--------------------------------
dirout = dirin
if not os.path.exists(dirout):
    os.makedirs(dirout)
verbose = 'T'
ACF = ['F','F']

t_burnin = []
t_end = []
## BB
bb_log = pathlib.Path(dirout+'parlog_fit_bb'+h5ext)
if bb_log.exists():
    Nmcmc = read_hdf5(dirout+'parlog_fit_bb', 'Last index')[0]
else:
    Nmcmc = read_hdf5(dirin+'input_master', 'Nmcmc')[0]
t_burnin.append(int( Nmcmc * 0.2 ))
t_end.append( Nmcmc )
## HB
hb_log = pathlib.Path(dirout+'parlog_fit_hb'+h5ext)
if hb_log.exists():
    Nmcmc = read_hdf5(dirout+'parlog_fit_hb', 'Last index')[0]
else:
    Nmcmc = read_hdf5(dirin+'input_master', 'Nmcmc')[0]
t_burnin.append(int( Nmcmc * 0.2 ))
t_end.append( Nmcmc )

## Write HDF5
##------------
write_hdf5(h5_analysis, 'program', [program], verbose=True)
write_hdf5(h5_analysis, 'output dir', [dirout], append=True, verbose=noisy)
write_hdf5(h5_analysis, 'verbose', [verbose], append=True, verbose=noisy)
write_hdf5(h5_analysis, 'ACF', ACF, append=True, verbose=noisy)
write_hdf5(h5_analysis, 't_burnin', t_burnin, append=True, verbose=noisy)
write_hdf5(h5_analysis, 't_end', t_end, append=True, verbose=noisy)
