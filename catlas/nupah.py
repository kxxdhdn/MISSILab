#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Nearby Universe PAH catalog

"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io import fits
import healpy as hp

from myfunclib import closest

cm = mpl.cm.gist_heat
# norm = mpl.colors.Normalize(vmin=0, vmax=.2)
norm = mpl.colors.SymLogNorm(linthresh=.01, vmin=0, vmax=10.)

Nside = 8192
Npix = hp.nside2npix(Nside) # Npix = 12 * Nside**2
res = hp.nside2resol(Nside, arcmin=True) * 60.
print('Nside = {} \n Npix = {} \n Resol = {} arcsec'.format(
	Nside, Npix, res))

fig = plt.figure(1, figsize=(24,10))

## WISE 12um all sky map
##-----------------------
mpath = '/Users/dhu/Downloads/WISE3/'
mfile = mpath+'wssa_sample_'+str(Nside)+'-bintable.fits'

hdul = fits.open(mfile)
all_sky_map = hdul[1].data.I12

## sample
##--------
spath = '/Users/dhu/Github/MISSILab/catlas/'
sfile = spath+'slist.txt'

dat = []
with open(sfile, 'r') as f:
	for line in f:
		dat.append(line.split())
data = Table(rows=dat, names=('TG_NAME', 'RA', 'DEC'), \
			 dtype=('S15', 'f8', 'f8'), masked=True)

## mollweide projection
##----------------------
hp.mollview(all_sky_map, fig=1, cmap=cm, 
	norm='symlog', 
	min=-5.e0, 
	max=5.e1, 
	bgcolor='lightgray',
	# xsize=2,
	# hold=True, 
	# return_projected_map=True
	)
hp.projscatter(data['RA'], data['DEC'], lonlat=True, 
	marker='*', c='y', 
	s=80.,
	)

# fig.savefig(spath+'skymap.png', transparent=True)

plt.show()
