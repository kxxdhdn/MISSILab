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

import os
import math
import numpy as np
import subprocess as SP

## laputan
from laputan.arrays import closest
from laputan.inout import read_fits, write_hdf5, read_hdf5

## local
from librarian import (croot, mroot,
                       res, TABLine, TABand, partuning)
         
## Path
##------
h5_path = mroot+'out/set_input'
dirin = croot+'../out/' # define input dir
if not os.path.exists(dirin):
    os.makedirs(dirin)
write_hdf5(h5_path, 'input dir', [dirin], verbose=True)

h5_analysis = dirin+'input_analysis'

program = 'analysis'
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
input_fit = croot+'input_sim_hb.py'
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
ACF = 'T'

Nmcmc = read_hdf5(dirin+'input_master', 'Nmcmc')[0]
t_burnin = int(Nmcmc*.1)
t_end = Nmcmc

## Write HDF5
##------------
write_hdf5(h5_analysis, 'program', [program], verbose=True)
write_hdf5(h5_analysis, 'output dir', [dirout], append=True, verbose=noisy)
write_hdf5(h5_analysis, 'verbose', [verbose], append=True, verbose=noisy)
write_hdf5(h5_analysis, 'ACF', [ACF], append=True, verbose=noisy)
write_hdf5(h5_analysis, 't_burnin', [t_burnin], append=True, verbose=noisy)
write_hdf5(h5_analysis, 't_end', [t_end], append=True, verbose=noisy)
