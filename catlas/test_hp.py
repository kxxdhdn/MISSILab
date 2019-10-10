#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Nearby Universe PAH catalog (healpy version)

"""
import numpy as np
import healpy as hp
import matplotlib as mpl
import matplotlib.pyplot as plt
# from astropy.table import Table
# import astropy.coordinates as coord
# import astropy.units as u
# from astropy.io import fits

# cmap = mpl.cm.viridis
cmap = mpl.cm.gist_heat
norm = mpl.colors.Normalize(vmin=0, vmax=1)
# norm = mpl.colors.LogNorm(vmin=1.e-2, vmax=1)
# norm = mpl.colors.SymLogNorm(thresh=1.e-2, vmin=0, vmax=1)

Nside = 4096
path = '/Users/dhu/Downloads/'
file = path+'wssa_sample_'+str(Nside)+'-bintable.fits'

Npix = hp.nside2npix(Nside)

all_sky_map = hp.read_map(file)

hp.mollview(all_sky_map,
	# coord=["G", "E"],
	# title="",
	cmap=cmap,
	norm="symlog",
	# min=1.e-2,
	# max=1,
	)
hp.graticule() # add meridians and parallels

plt.show()
