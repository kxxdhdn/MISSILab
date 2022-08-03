#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Unit library

"""

import numpy as np

from iolib import read_fits, write_fits
from astrolib import pix2sr, get_pc

def Jy_per_pix_to_MJy_per_sr(filIN, filOUT=None):
    '''
    Convert image unit from Jy/pix to MJy/sr
    '''
    hd = read_fits(filIN)
    oldimage = hd.data

    else:
    cdelt = get_pc(header=hd.header).cdelt
    ## gmean( Jy/MJy / sr/pix )
    ufactor = np.sqrt(np.prod(1.e-6/pix2sr(1., cdelt)))
    # print(cdelt, unit_fac)
    newimage = oldimage * ufactor
    hdr = swp.refheader
    hdr['BUNIT'] = 'MJy/sr'
    write_fits(path_cal+src+'_'+phot, swp.refheader, image_phot)

    return newimage
