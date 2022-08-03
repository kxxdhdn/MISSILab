#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Available variables and functions:

  PATH: croot, mroot
  DATA: res, TABLine, TABand, WAVLine, WAVBand
  FUNC: plotname, plotcorr, calcorr
  FUNC: partuning, calexpansion, calextraction

"""

import os
import numpy as np


## Path of current file
croot = os.path.dirname(os.path.abspath(__file__))+'/'
## MILES root path
mroot = '/Users/dhu/ownCloud/MILES/'


## Parameter names for plot
def plotname(parname):
    '''
    ------ INPUT ------
    parname             Parameter name in the model
    ------ OUTPUT ------
    plotname            Parameter name in the plot
    '''
    if parname=='lnAv1':
        plotname = r'ln($A_V$)'
    elif parname=='lnLstar1':
        plotname = r'ln($L^{\star}$)'
    elif parname=='lnFcont1':
        plotname = r'ln($F^{BE\ Z96}_{15\ \mu m}$)'
    elif parname=='lnT1':
        plotname = r'ln($T^{BE\ Z96}$)'
    elif parname=='lnFcont2':
        plotname = r'ln($F^{Gra\ D03}_{15\ \mu m}$)'
    elif parname=='lnT2':
        plotname = r'ln($T^{Gra\ D03}$)'
    elif parname=='lnFcont3':
        plotname = r'ln($F^{Sil\ D03}_{15\ \mu m}$)'
    elif parname=='lnT3':
        plotname = r'ln($T^{Sil\ D03}$)'
    elif parname=='lnRline1':
        plotname = r'ln($R^{Bra}$)'
    elif parname=='lnRline2':
        plotname = r'ln($R^{ArII}$)'
    elif parname=='lnRline3':
        plotname = r'ln($R^{ArIII1}$)'
    elif parname=='lnRline4':
        plotname = r'ln($R^{H2S3}$)'
    elif parname=='lnRline5':
        plotname = r'ln($R^{SIV}$)'
    elif parname=='lnRline6':
        plotname = r'ln($R^{H2S2}$)'
    elif parname=='lnRline7':
        plotname = r'ln($R^{NeII}$)'
    elif parname=='lnRline8':
        plotname = r'ln($R^{NeIII1}$)'
    elif parname=='lnRline9':
        plotname = r'ln($R^{H2S1}$)'
    elif parname=='lnRline10':
        plotname = r'ln($R^{SIII1}$)'
    elif parname=='lnRband1':
        plotname = r'ln($R^{3.3}$)'
    elif parname=='lnRband2':
        plotname = r'ln($R^{3.4}$)'
    elif parname=='lnRband3':
        plotname = r'ln($R^{3.5}$)'
    elif parname=='lnRband4':
        plotname = r'ln($R^{6.2\ (1)}$)'
    elif parname=='lnRband5':
        plotname = r'ln($R^{6.2\ (2)}$)'
    elif parname=='lnRband6':
        plotname = r'ln($R^{7.7\ (p)}$)'
    elif parname=='lnRband7':
        plotname = r'ln($R^{7.7\ (1)}$)'
    elif parname=='lnRband8':
        plotname = r'ln($R^{7.7\ (2)}$)'
    elif parname=='lnRband9':
        plotname = r'ln($R^{8.6}$)'
    elif parname=='lnRband10':
        plotname = r'ln($R^{11.0}$)'
    elif parname=='lnRband11':
        plotname = r'ln($I^{11.2}$)'
    elif parname=='lnRband12':
        plotname = r'ln($R^{11.3\ (p)}$)'
    elif parname=='lnRband13':
        plotname = r'ln($R^{12.7\ (1)}$)'
    elif parname=='lnRband14':
        plotname = r'ln($R^{12.7\ (2)}$)'
    elif parname=='lnRband15':
        plotname = r'ln($R^{17.0\ (p)}$)'
    else:
        plotname = parname

    return plotname

## Parameter names for plot
def plotcorr(corrname):
    '''
    ------ INPUT ------
    corrname             Correlation name in the model
    ------ OUTPUT ------
    plotcorr             Correlation name in the plot
    '''
    if corrname=='I3.3_ov_I11.3':
        plotcorr = 'I3.3/I11.3'
    elif corrname=='I3.4_ov_I11.3':
        plotcorr = 'I3.4/I11.3'
    elif corrname=='I6.2_ov_I11.3':
        plotcorr = 'I6.2/I11.3'
    elif corrname=='I7.7_ov_I11.3':
        plotcorr = 'I7.7/I11.3'
    elif corrname=='I12.7_ov_I11.3':
        plotcorr = 'I12.7/I11.3'

    return plotcorr

## calculate
def calcorr(pararr, corrname, indpar, parerr=None, chi2=False):
    '''
    ------ INPUT ------
    pararr               Parameter array
    corrname             Correlation name in the model
    indpar               Default index list of par used for corr
    chi2                 If true, calculated propagated errors
    ------ OUTPUT ------
    calcorr              Index list of par used for corr
    '''
    if corrname=='I3.3_ov_I11.3':
        indcorr = [1]
    elif corrname=='I3.4_ov_I11.3':
        indcorr = [2,3]
    elif corrname=='I6.2_ov_I11.3':
        indcorr = [4,5]
    elif corrname=='I7.7_ov_I11.3':
        indcorr = [6,7,8]
    elif corrname=='I12.7_ov_I11.3':
        indcorr = [9,10]

    v0 = 1+np.exp(pararr[indpar[0]])
    
    calcorr = 0
    for i in indcorr:
        calcorr += np.exp(pararr[indpar[i]]) / v0
        
    if chi2:
        e0 = ( np.exp(pararr[indpar[0]])*parerr[indpar[0]] )**2
        val = 0
        val2 = 0
        for i in indcorr:
            val += ( np.exp(pararr[indpar[i]])*parerr[indpar[i]] )**2
            val2 += np.exp(pararr[indpar[i]])
        correrr = np.sqrt(val + val2**2 * e0) / v0
    else:
        correrr = None

    return calcorr, correrr

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
                   ('label','FeVII1'),('wave',7.8145) ]),
            dict([ ('name','[Ar!E !NV] 3P!U1!N-3P!U2!N'),
                   ('label','ArV1'),('wave',7.90158) ]),
            dict([ ('name','H!D2!N 0-0 S(4)'),
                   ('label','H2S4'),('wave',8.02505) ]),
            dict([ ('name','H!E !NI 7-10'),
                   ('label','HI7-10'),('wave',8.760064) ]),
            dict([ ('name','[Ar!E !NIII] 3P!D2!N-3P!D1!N'),
                   ('label','ArIII1'),('wave',8.99103) ]),
            dict([ ('name','[Ni!E !NVI] 4P!U5/2!N-4F!U5/2!N'),
                   ('label','NiVI'),('wave',9.042) ]), 
            dict([ ('name','[Fe!E !NVII] 3F!D2!N-3F!D3!N'),
                   ('label','FeVII2'),('wave',9.5267) ]),
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
                   ('label','ArV2'),('wave',13.10219) ]),
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

WAVLine = [ 4.052    , # Bra
            5.128657 , # HI6-10
            5.51116  , # H2S7
            5.6098   , # MgV1
            5.908213 , # Huc
            5.981    , # KIV
            6.10856  , # H2S6
            6.709    , # ClV
            6.90952  , # H2S5
            6.947984 , # HeII1
            6.985274 , # ArII
            7.3178   , # NaIII
            7.459858 , # Pfa
            7.502493 , # Hub
            7.6524   , # NeVI
            7.8145   , # FeVII1
            7.90158  , # ArV1
            8.02505  , # H2S4
            8.760064 , # HI7-10
            8.99103  , # ArIII1
            9.042    , # NiVI
            9.5267   , # FeVII2
            9.66491  , # H2S3
            9.713475 , # HeII2
            10.3385  , # SiI
            10.5105  , # SIV
            12.27861 , # H2S2
            12.368527, # Hua
            12.81355 , # NeII
            13.10219 , # ArV2
            13.521   , # MgV2
            14.32168 , # NeV1
            15.555   , # NeIII1
            17.03483 , # H2S1
            17.608246, # HI11-18
            18.7129  , # SIII1
            19.061898, # HI7-8
            21.8291  , # ArIII2
            24.3175  , # NeV2
            25.8903  , # OIV
            25.9882  , # FeII1
            28.21883 , # H2S0
            33.4810  , # SIII2
            34.8152  , # SiII
            35.3491  , # FeII2
            36.0135  , ] # NeIII2

WAVBand = [ 3.291    ,
            3.399    ,
            3.499    ,
            5.2394667,
            5.6437461,
            5.7490305,
            6.0106598, # Small 6.0
            6.2034273,
            6.2672596,
            6.6273788,
            6.8548833,
            7.0791725, # Small 7.1
            7.6000000,
            7.6171123,
            7.8704769,
            8.3623706, # Small 8.3
            8.6204540,
            9.5244838,
            10.707132,
            11.038349,
            11.237893,
            11.400432, # Plateau 11.3
            11.796389,
            11.949674,
            12.626842,
            12.760273,
            13.559342,
            14.257133,
            15.893117,
            16.482868,
            17.082868,
            17.428485,
            17.771096,
            18.925630, ]

def partuning(dictune, labQ, labL, labB, labE, labS):
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
    ds                  output properties
      name                parinfo name (list)
      comp                parinfo comp
      fixed               parinfo fixed
      limited             parinfo limited
      limits              parinfo limits
      model               parinfo model
      hyper               parinfo hyper
      tied                parinfo tied
      value               parinfo value
    ------ OUTPUT ------
    '''
    
    Ncont = len(labQ)
    Nline = len(labL)
    Nband = len(labB)
    Nextc = len(labE)
    Nstar = len(labS)
    Npar = 2*Ncont + 3*Nline + 4*Nband + Nextc + Nstar
    
    ## Set indL (index for labL)
    indL = []
    for lab in labL:
        for tabL in TABLine:
            if tabL['label']==lab.rstrip():
                indL.append(TABLine.index(tabL))
    ## Set indB (index for labB)
    indB = []
    for lab in labB:
        for tabB in TABand:
            if tabB['label']==lab.rstrip():
                indB.append(TABand.index(tabB))

    ## Default setting
    name = np.empty((Npar,), dtype=('<U30')) # str length <30
    namall = np.empty((Npar,), dtype=('<U30')) # specific for partuning
    comp = np.empty((Npar,), dtype=('<U30'))
    fixed = np.array(['T' for i in range(Npar)]) # fixed by default
    limited = np.array([('F','F') for i in range(Npar)])
    limits = np.array([(0.,0.) for i in range(Npar)])
    model = np.array(['T' for i in range(Npar)])
    hyper = np.array(['F' for i in range(Npar)])
    tied = np.empty((Npar,), dtype=('<U30'))
    value = np.array([0. for i in range(Npar)])

    ## Parinfo assignment
    ##--------------------
    indict = {'lnFcont':[],'lnT':[],
              'lnRline':[],'Cline':[],'Wline':[],
              'lnRband':[],'Cband':[],'WSband':[],'WLband':[],
              'lnAv':[],
              'lnLstar':[],} # ordered parameter number list for namall cases
    i0 = 0
    ## Extc
    for i in range(Nextc):
        indict['lnAv'].append(i0+i)
        # if len(str(i+1))==1:
        #     name[i0+i] = 'lnAv0'+str(i+1)
        # else:
        name[i0+i] = 'lnAv'+str(i+1)
        namall[i0+i] = 'lnAv'
        value[i0+i] = 0. # 1 [mag]
        comp[i0+i] = 'EXTC'
    i0 += Nextc
    ## Star
    for i in range(Nstar):
        indict['lnLstar'].append(i0+i)
        # if len(str(i+1))==1:
        #     name[i0+i] = 'lnLstar0'+str(i+1)
        # else:
        name[i0+i] = 'lnLstar'+str(i+1)
        namall[i0+i] = 'lnLstar'
        value[i0+i] = 0. # 1 [W/m2]
        comp[i0+i] = 'STAR'
    i0 += Nstar
    ## Cont
    for i in range(Ncont):
        indict['lnFcont'].append(i0+2*i)
        # if len(str(i+1))==1:
        #     name[i0+2*i] = 'lnFcont0'+str(i+1)
        # else:
        name[i0+2*i] = 'lnFcont'+str(i+1)
        namall[i0+2*i] = 'lnFcont'
        value[i0+2*i] = 0.
        indict['lnT'].append(i0+2*i+1)
        # if len(str(i+1))==1:
        #     name[i0+2*i+1] = 'lnT0'+str(i+1)
        # else:
        name[i0+2*i+1] = 'lnT'+str(i+1)
        namall[i0+2*i+1] = 'lnT'
        value[i0+2*i+1] = 4. # 54.60 [K]
    for i in range(2*Ncont):
        comp[i0+i] = 'CONT'
    i0 += 2*Ncont
    ## Line
    for i in range(Nline):
        indict['lnRline'].append(i0+3*i)
        # if len(str(i+1))==1:
        #     name[i0+3*i] = 'lnRline0'+str(i+1)
        # else:
        name[i0+3*i] = 'lnRline'+str(i+1)
        namall[i0+3*i] = 'lnRline'
        value[i0+3*i] = 0.
        indict['Cline'].append(i0+3*i+1)
        # if len(str(i+1))==1:
        #     name[i0+3*i+1] = 'Cline0'+str(i+1)
        # else:
        name[i0+3*i+1] = 'Cline'+str(i+1)
        namall[i0+3*i+1] = 'Cline'
        value[i0+3*i+1] = TABLine[indL[i]]['wave']
        indict['Wline'].append(i0+3*i+2)
        # if len(str(i+1))==1:
        #     name[i0+3*i+2] = 'Wline0'+str(i+1)
        # else:
        name[i0+3*i+2] = 'Wline'+str(i+1)
        namall[i0+3*i+2] = 'Wline'
        for r in res:
            if value[i0+3*i+1]<5.:
                if r['name']=='AKARI_NG':
                    value[i0+3*i+2] = r['dwovw'] * value[i0+3*i+1]
            elif value[i0+3*i+1]<14.29:
                if r['name']=='SL':
                    value[i0+3*i+2] = r['dwovw'] * value[i0+3*i+1]
            else:
                if r['name']=='LL':
                    value[i0+3*i+2] = r['dwovw'] * value[i0+3*i+1]
    for i in range(3*Nline):
        comp[i0+i] = 'LINE'
    i0 += 3*Nline
    ## Band
    for i in range(Nband):
        indict['lnRband'].append(i0+4*i)
        # if len(str(i+1))==1:
        #     name[i0+4*i] = 'lnRband0'+str(i+1)
        # else:
        name[i0+4*i] = 'lnRband'+str(i+1)
        namall[i0+4*i] = 'lnRband'
        value[i0+4*i] = 0.
        indict['Cband'].append(i0+4*i+1)
        # if len(str(i+1))==1:
        #     name[i0+4*i+1] = 'Cband0'+str(i+1)
        # else:
        name[i0+4*i+1] = 'Cband'+str(i+1)
        namall[i0+4*i+1] = 'Cband'
        value[i0+4*i+1] = TABand[indB[i]]['wave']
        indict['WSband'].append(i0+4*i+2)
        # if len(str(i+1))==1:
        #     name[i0+4*i+2] = 'WSband0'+str(i+1)
        # else:
        name[i0+4*i+2] = 'WSband'+str(i+1)
        namall[i0+4*i+2] = 'WSband'
        value[i0+4*i+2] = TABand[indB[i]]['sigmaS']
        indict['WLband'].append(i0+4*i+3)
        # if len(str(i+1))==1:
        #     name[i0+4*i+3] = 'WLband0'+str(i+1)
        # else:
        name[i0+4*i+3] = 'WLband'+str(i+1)
        namall[i0+4*i+3] = 'WLband'
        value[i0+4*i+3] = TABand[indB[i]]['sigmaL']
        # for r in res:
        #     if value[i0+4*i+1]<5.:
        #         if r['name']=='AKARI_NG':
        #             width = r['dwovw'] * value[i0+4*i+1]
        #             if value[i0+4*i+2]<width:
        #                 value[i0+4*i+2] = width
        #             if value[i0+4*i+3]<width:
        #                 value[i0+4*i+3] = width
        #     elif value[i0+4*i+1]<14.29:
        #         if r['name']=='SL':
        #             width = r['dwovw'] * value[i0+4*i+1]
        #             if value[i0+4*i+2]<width:
        #                 value[i0+4*i+2] = width
        #             if value[i0+4*i+3]<width:
        #                 value[i0+4*i+3] = width
        #     else:
        #         if r['name']=='LL':
        #             width = r['dwovw'] * value[i0+4*i+1]
        #             if value[i0+4*i+2]<width:
        #                 value[i0+4*i+2] = width
        #             if value[i0+4*i+3]<width:
        #                 value[i0+4*i+3] = width
    for i in range(4*Nband):
        comp[i0+i] = 'BAND'
    i0 += 4*Nband
    
    ## Param tuning
    ##--------------
    for tune in dictune:
        ind = []
        if 'namall' in tune:
            if tune['namall']!='default':
                ind = indict[tune['namall']]
                
        if 'name' in tune:
            n = tune['name']
            # n = tune['name'].encode('ascii','ignore')
            if n in name:
                ind = [ list(name).index(n) ]

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
                value[i] = tune['value'] # Personalized value

    ## Outputs
    ds = type('', (), {})()
    ds.name = name
    ds.comp = comp
    ds.fixed = fixed
    ds.limited = limited
    ds.limits = limits
    ds.model = model
    ds.hyper = hyper
    ds.tied = tied
    ds.value = value

    return ds

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
