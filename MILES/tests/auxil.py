#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Available variables and functions:

  PATH: croot, mroot
  DATA: res, TABLine, TABand
  FUNC: partuning, calexpansion, calextraction

"""

import os
import numpy as np


## Path of current file
croot = os.path.dirname(os.path.abspath(__file__))+'/'
## MILES root path
mroot = '/Users/dhu/Github/MISSILE/MILES/'


res = [ dict([ ('name','CAM'),('dwovw',0.010373835) ]),
        dict([ ('name','SL'),('dwovw',0.0055722841) ]),
        dict([ ('name','SH'),('dwovw',0.00075127093) ]),
        dict([ ('name','LL'),('dwovw',0.0057517091) ]),
        dict([ ('name','LH'),('dwovw',0.00070906159) ]),
        dict([ ('name','SWS'), ('dwovw',0.0011044189) ]),
        dict([ ('name','SWSfine'),('dwovw',0.00036671469) ]),
        dict([ ('name','AKARI_NG'),('dwovw',0.00881818) ]) ]

TABLine = [ dict([ ('name','Bracket alpha'),
                   ('label','Bra'),('wave',4.052) ]),
            dict([ ('name','H!E !NI 6-10'),
                   ('label','HI6-10'),('wave',5.128657) ]), # AKARI_NG (2)
            dict([ ('name','H!D2!N 0-0 S(7)'),
                   ('label','H2S7'),('wave',5.51116) ]),
            dict([ ('name','[Mg!E !NV] 3P!U2!N-3P!U1!N'),
                   ('label','MgV1'),('wave',5.6098) ]),
            dict([ ('name','Humphreys gamma'),
                   ('label','Huc'),('wave',5.908213) ]),
            dict([ ('name','[K!E !NIV]'),
                   ('label','KIV'),('wave',5.981) ]), 
            dict([ ('name','H!D2!N 0-0 S(6)'),
                   ('label','H2S6'),('wave',6.10856) ]), 
            dict([ ('name','[Cl!E !NV] 2P!U0!N-2P!U0!N'),
                   ('label','ClV'),('wave',6.709) ]),
            dict([ ('name','H!D2!N 0-0 S(5)'),
                   ('label','H2S5'),('wave',6.90952) ]), 
            dict([ ('name','HeII 8-9'),
                   ('label','HeII1'),('wave',6.947984) ]),
            dict([ ('name','[Ar!E !NII] 2P!D3/2!N-2P!D1/2!N'),
                   ('label','ArII'),('wave',6.985274) ]),
            dict([ ('name','[Na!E !NIII] 2P!U0!N-2P!U0!N'),
                   ('label','NaIII'),('wave',7.3178) ]),
            dict([ ('name','Pfund alpha'),
                   ('label','Pfa'),('wave',7.459858) ]),
            dict([ ('name','Humphreys beta'),
                   ('label','Hub'),('wave',7.502493) ]), # SL2 (12)
            dict([ ('name','[Ne!E !NVI] 2P!U0!N-2P!U0!N'),
                   ('label','NeVI'),('wave',7.6524) ]),
            dict([ ('name','[Fe!E !NVII] 3F!U3!N-3F!U4!N'),
                   ('label','FeVII'),('wave',7.8145) ]),
            dict([ ('name','[Ar!E !NV] 3P!U1!N-3P!U2!N'),
                   ('label','ArV'),('wave',7.90158) ]),
            dict([ ('name','H!D2!N 0-0 S(4)'),
                   ('label','H2S4'),('wave',8.02505) ]),
            dict([ ('name','H!E !NI 7-10'),
                   ('label','HI7-10'),('wave',8.760064) ]),
            dict([ ('name','[Ar!E !NIII] 3P!D2!N-3P!D1!N'),
                   ('label','ArIII1'),('wave',8.99103) ]),
            dict([ ('name','[Ni!E !NVI] 4P!U5/2!N-4F!U5/2!N'),
                   ('label','NiVI'),('wave',9.042) ]), 
            dict([ ('name','[Fe!E !NVII] 3F!D2!N-3F!D3!N'),
                   ('label','FeVII'),('wave',9.5267) ]),
            dict([ ('name','H!D2!N 0-0 S(3)'),
                   ('label','H2S3'),('wave',9.66491) ]),
            dict([ ('name','He!E !NII 9-10'),
                   ('label','HeII2'),('wave',9.713475) ]),
            dict([ ('name','[Si!E !NI] 1P!D1!N-1P!D2!N'),
                   ('label','SiI'),('wave',10.3385) ]),
            dict([ ('name','[S!E !NIV] 2P!D3/2!N-2P!D1/2!N'),
                   ('label','SIV'),('wave',10.5105) ]),
            dict([ ('name','H!D2!N 0-0 S(2)'),
                   ('label','H2S2'),('wave',12.27861) ]),
            dict([ ('name','Humphreys alpha'),
                   ('label','Hua'),('wave',12.368527) ]),
            dict([ ('name','[Ne!E !NII] 2P!D3/2!N-2P!D1/2!N'),
                   ('label','NeII'),('wave',12.81355) ]),
            dict([ ('name','[Ar!E !NV] 3P!D0!N-3P!D1!N'),
                   ('label','ArV'),('wave',13.10219) ]),
            dict([ ('name','[Mg!E !NV] 3P!D1!N-3P!D0!N'),
                   ('label','MgV2'),('wave',13.521) ]), # SL1 (17)
            dict([ ('name','[Ne!E !NV] 3P!D1!N-3P!D2!N'),
                   ('label','NeV1'),('wave',14.32168) ]),
            dict([ ('name','[Ne!E !NIII] 3P!D2!N-3P!D1!N'),
                   ('label','NeIII1'),('wave',15.555) ]),
            dict([ ('name','H!D2!N 0-0 S(1)'),
                   ('label','H2S1'),('wave',17.03483) ]),
            dict([ ('name','H!E !NI 11-18'),
                   ('label','HI11-18'),('wave',17.608246) ]),
            dict([ ('name','[S!E !NIII] 3P!D2!N-3P!D1!N'),
                   ('label','SIII1'),('wave',18.7129) ]),
            dict([ ('name','H!E !NI 7-8'),
                   ('label','HI7-8'),('wave',19.061898) ]),
            dict([ ('name','[Ar!E !NIII] 3P!U1!N-3P!U0!N'),
                   ('label','ArIII2'),('wave',21.8291) ]), # LL2 (6)
            dict([ ('name','[Ne!E !NV] 3P!U0!N-3P!U1!N'),
                   ('label','NeV2'),('wave',24.3175) ]),
            dict([ ('name','[O!E !NIV] 2P!D3/2!N-2P!D1/2!N'),
                   ('label','OIV'),('wave',25.8903) ]),
            dict([ ('name','[Fe!E !NII] 6D!D7/2!N-6D!D9/2!N'),
                   ('label','FeII1'),('wave',25.9882) ]),
            dict([ ('name','H!D2!N 0-0 S(0)'),
                   ('label','H2S0'),('wave',28.21883) ]),
            dict([ ('name','[S!E !NIII] 3P!D1!N-3P!D0!N'),
                   ('label','SIII2'),('wave',33.4810) ]),
            dict([ ('name','[Si!E !NII] 2P!D3/2!N-2P!D1/2!N'),
                   ('label','SiII'),('wave',34.8152) ]),
            dict([ ('name','[Fe!E !NII] 6D!D5/2!N-6D!D7/2!N'),
                   ('label','FeII2'),('wave',35.3491) ]),
            dict([ ('name','[Ne!E !NIII] 3P!D1!N-3P!D0!N'),
                   ('label','NeIII2'),('wave',36.0135) ]) ] # LL1 (9)

TABand = [ dict([ ('label','Main 3.3'),('wave',3.291),
                  ('sigmaS',0.020),('sigmaL',0.019) ]),
           dict([ ('label','Main 3.4'),('wave',3.399),
                  ('sigmaS',0.011),('sigmaL',0.024) ]),
           dict([ ('label','Small 3.5'),('wave',3.499),
                  ('sigmaS',0.077),('sigmaL',0.071) ]), # AKARI_NG (3)
           dict([ ('label','Small 5.2'),('wave',5.2394667),
                  ('sigmaS',0.025218240),('sigmaL',0.058333420) ]),
           dict([ ('label','Small 5.7 (1)'),('wave',5.6437461),
                  ('sigmaS',0.040000000),('sigmaL',0.080000000) ]),
           dict([ ('label','Small 5.7 (2)'),('wave',5.7490305),
                  ('sigmaS',0.040000000),('sigmaL',0.080000000) ]),
           dict([ ('label','Small 6.0'),('wave',6.0106598),
                  ('sigmaS',0.040000000),('sigmaL',0.066646401) ]),
           dict([ ('label','Main 6.2 (1)'),('wave',6.2034273),
                  ('sigmaS',0.031300317),('sigmaL',0.060000000) ]),
           dict([ ('label','Main 6.2 (2)'),('wave',6.2672596),
                  ('sigmaS',0.036922155),('sigmaL',0.11633640) ]),
           dict([ ('label','Small 6.6'),('wave',6.6273788),
                  ('sigmaS',0.12000000),('sigmaL',0.12000000) ]),
           dict([ ('label','Small 6.8'),('wave',6.8548833),
                  ('sigmaS',0.080000000),('sigmaL',0.080000000) ]),
           dict([ ('label','Small 7.1'),('wave',7.0791725),
                  ('sigmaS',0.080000000),('sigmaL',0.080000000) ]), # SL2 (9)
           dict([ ('label','Plateau 7.7'),('wave',7.6000000),
                  ('sigmaS',0.48000000),('sigmaL',0.50247515) ]),
           dict([ ('label','Main 7.7 (1)'),('wave',7.6171123),
                  ('sigmaS',0.11856752),('sigmaL',0.14531480) ]),
           dict([ ('label','Main 7.7 (2)'),('wave',7.8704769),
                  ('sigmaS',0.16998357),('sigmaL',0.24523967) ]),
           dict([ ('label','Small 8.3'),('wave',8.3623706),
                  ('sigmaS',0.016256724),('sigmaL',0.016256724) ]),
           dict([ ('label','Main 8.6'),('wave',8.6204540),
                  ('sigmaS',0.18340577),('sigmaL',0.13337054) ]),
           dict([ ('label','Small 9.5'),('wave',9.5244838),
                  ('sigmaS',0.10766965),('sigmaL',0.60000000) ]),
           dict([ ('label','Small 10.7'),('wave',10.707132),
                  ('sigmaS',0.10000000),('sigmaL',0.10000000) ]),
           dict([ ('label','Small 11.0'),('wave',11.038349),
                  ('sigmaS',0.026989462),('sigmaL',0.073146141) ]),
           dict([ ('label','Main 11.2'),('wave',11.237893),
                  ('sigmaS',0.053507232),('sigmaL',0.15254629) ]),
           dict([ ('label','Plateau 11.3'),('wave',11.400432),
                  ('sigmaS',0.72000000),('sigmaL',0.63657944) ]),
           dict([ ('label','Small 11.8'),('wave',11.796389),
                  ('sigmaS',0.020813349),('sigmaL',0.020813349) ]),
           dict([ ('label','Small 11.9'),('wave',11.949674),
                  ('sigmaS',0.080352283),('sigmaL',0.22192473) ]),
           dict([ ('label','Main 12.7 (1)'),('wave',12.626842),
                  ('sigmaS',0.20000000),('sigmaL',0.094424464) ]),
           dict([ ('label','Main 12.7 (2)'),('wave',12.760273),
                  ('sigmaS',0.080436118),('sigmaL',0.14000000) ]),
           dict([ ('label','Small 13.6'),('wave',13.559342),
                  ('sigmaS',0.15954880),('sigmaL',0.16054015) ]),
           dict([ ('label','Small 14.2'),('wave',14.257133),
                  ('sigmaS',0.15208135),('sigmaL',0.058951597) ]), # SL1 (16)
           dict([ ('label','Small 15.6'),('wave',15.893117),
                  ('sigmaS',0.17857214),('sigmaL',0.20000000) ]),
           dict([ ('label','Small 16.4'),('wave',16.482868),
                  ('sigmaS',0.10000000),('sigmaL',0.058462024) ]),
           dict([ ('label','Plateau 17.0'),('wave',17.082868),
                  ('sigmaS',0.49775906),('sigmaL',0.56119197) ]),
           dict([ ('label','Small 17.4'),('wave',17.428485),
                  ('sigmaS',0.10000000),('sigmaL',0.10000000) ]),
           dict([ ('label','Small 17.8'),('wave',17.771096),
                  ('sigmaS',0.030799172),('sigmaL',0.075249330) ]),
           dict([ ('label','Small 18.9'),('wave',18.925630),
                  ('sigmaS',0.034553879),('sigmaL',0.11570587) ]) ] # LL2 (6)

def partuning(dictune, Ncont, Nline, Nband, Nextc,
              name, fixed, limited, limits, model, hyper, tied, value):
    '''
    ------ INPUT ------
    dictune             dict of param tuning info
      name                name of param (covering 'namall' if co-exist)
      namall              tune all param of the same type
      fixed               'T'/'F'
      limited             ('T'/'F','T'/'F')
      limits              (inf,sup) float tuple
      hyper               'T'/'F'
      tied                name of tied param
      value               float type
    Ncont               N of continuum components
    Nline               N of lines
    Nband               N of bands
    Nextc               N of extinction compo
    name                parinfo name (list)
    fixed               parinfo fixed
    limited             parinfo limited
    limits              parinfo limits
    model               parinfo model
    hyper               parinfo hyper
    tied                parinfo tied
    value               parinfo value
    ------ OUTPUT ------
    '''
    for tune in dictune:
        if 'namall' in tune:
            ind = []
            if tune['namall']=='lnFcont':
                i0 = 0
                ind = [i0+2*i for i in range(Ncont)]
            elif tune['namall']=='lnT':
                i0 = 0
                ind = [i0+2*i+1 for i in range(Ncont)]
            elif tune['namall']=='lnRline':
                i0 = 2*Ncont
                ind = [i0+3*i for i in range(Nline)]
            elif tune['namall']=='Cline':
                i0 = 2*Ncont
                ind = [i0+3*i+1 for i in range(Nline)]
            elif tune['namall']=='Wline':
                i0 = 2*Ncont
                ind = [i0+3*i+2 for i in range(Nline)]
            elif tune['namall']=='lnRband':
                i0 = 2*Ncont + 3*Nline
                ind = [i0+4*i for i in range(Nband)]
            elif tune['namall']=='Cband':
                i0 = 2*Ncont + 3*Nline
                ind = [i0+4*i+1 for i in range(Nband)]
            elif tune['namall']=='WSband':
                i0 = 2*Ncont + 3*Nline
                ind = [i0+4*i+2 for i in range(Nband)]
            elif tune['namall']=='WLband':
                i0 = 2*Ncont + 3*Nline
                ind = [i0+4*i+3 for i in range(Nband)]
            elif tune['namall']=='lnAv':
                i0 = 2*Ncont + 3*Nline + 4*Nband
                ind = [i0+i for i in range(Nextc)]
            elif tune['namall']=='lnFstar':
                i0 = 2*Ncont + 3*Nline + 4*Nband + Nextc
                ind = [i0]
        if 'name' in tune:
            n = tune['name']
            # n = tune['name'].encode('ascii','ignore')
            if n in name:
                ind = [ list(name).index(n) ]
                if 'fixed' in tune:
                    fixed[i] = tune['fixed']
        
        for i in ind:
            if 'fixed' in tune:
                fixed[i] = tune['fixed']
            if 'limited' in tune:
                limited[i] = tune['limited']
            if 'limits' in tune:
                limits[i] = tune['limits']
            if 'model' in tune:
                model[i] = tune['model']
            if 'hyper' in tune:
                hyper[i] = tune['hyper']
            if 'tied' in tune:
                tied[i] = tune['tied']
            if 'value' in tune:
                value[i] = tune['value']

def calexpansion(calib, wvl, instr):
    '''
    Attribute calibration/scale factors to
    wavelength grid according to given instruments
    ------ INPUT ------
    calib               list of calibration factors
    x                   wavelength grid
    instr               list of instruments
    ------ OUTPUT ------
    carray              calib factor array (same size of x)
    '''
    carray = np.zeros(wvl.shape)
    for i, ins in enumerate(instr):
        if ins=='IRC_NG':
            lim = [2.5, 5.0]
        elif ins=='IRS_SL2':
            lim = [5.21, 7.56]
        elif ins=='IRS_SL1':
            lim = [7.57, 14.28]
        elif ins=='IRS_LL2':
            lim = [14.29, 20.66]
        elif ins=='IRS_LL1':
            lim = [20.67, 38.00]
        elif ins=='IRS_SH':
            lim = [10.00, 19.19]
        elif ins=='IRS_LH':
            lim = [19.20, 37.10]

        for iw, w in enumerate(wvl):
            if w>lim[0] and w<lim[1]:
                carray[iw] = calib[i]

    return carray

def calextraction(carray, wvl, instr):
    '''
    Extract calibration/scale factors of the instruments
    ------ INPUT ------
    carray              enbedded calibration factors
    x                   wavelength grid
    instr               list of instruments
    ------ OUTPUT ------
    calib               calib factors (same size of instr)
    '''
    for i, ins in enumerate(instr):
        if ins=='IRC_NG':
            lim = [2.5, 5.0]
        elif ins=='IRS_SL2':
            lim = [5.21, 7.56]
        elif ins=='IRS_SL1':
            lim = [7.57, 14.28]
        elif ins=='IRS_LL2':
            lim = [14.29, 20.66]
        elif ins=='IRS_LL1':
            lim = [20.67, 38.00]
        elif ins=='IRS_SH':
            lim = [10.00, 19.19]
        elif ins=='IRS_LH':
            lim = [19.20, 37.10]

        for iw, w in enumerate(wvl):
            if w>lim[0] and w<lim[1]:
                calib[i] = carray[iw]
                break

    return calib
