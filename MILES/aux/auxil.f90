!******************************************************************************
!*
!*                              MILES Auxiliary Data
!*
!******************************************************************************


MODULE auxil

  USE utilities, ONLY: DP
  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: Ncont_max = 31, Nline_max = 46, Nband_max = 34
  INTEGER, PARAMETER, PUBLIC :: Nextc_max = 10,  Nstar_max = 1
  REAL(DP), PARAMETER, PUBLIC :: Cband_sig = 0.1_DP ! band center var range [micron]
  
  !!-------------------------
  !! Instrumental resolution
  !!-------------------------
  TYPE spect_res
    !! dlam/lam = dnu/nu
    REAL(DP) :: dw_w_CAM = 0.010373835_DP
    REAL(DP) :: dw_w_SWS = 0.0011044189_DP
    REAL(DP) :: dw_w_SWSfine = 0.00036671469_DP
    REAL(DP) :: dw_w_SL = 0.0055722841_DP
    REAL(DP) :: dw_w_SH = 0.00075127093_DP
    REAL(DP) :: dw_w_LL = 0.0057517091_DP
    REAL(DP) :: dw_w_LH = 0.00070906159_DP
    REAL(DP) :: dw_w_AKARI_NG = 0.00881818_DP
  END TYPE spect_res
  PUBLIC :: spect_res
  
  TYPE(spect_res), SAVE, PUBLIC :: res

  !!-----------------
  !! Spectral ranges
  !!-----------------
  TYPE spect_ran
    !! [w_start, w_end]
    REAL(DP), DIMENSION(2) :: wlim_AKARI_NG = (/ 2.50, 5.00 /)
    REAL(DP), DIMENSION(2) :: wlim_SL2 = (/ 5.21, 7.56 /)
    REAL(DP), DIMENSION(2) :: wlim_SL1 = (/ 7.57, 14.28 /)
    REAL(DP), DIMENSION(2) :: wlim_LL2 = (/ 14.29, 20.66 /)
    REAL(DP), DIMENSION(2) :: wlim_LL1 = (/ 20.67, 38.00 /)
    REAL(DP), DIMENSION(2) :: wlim_SH = (/ 10.00, 19.19 /)
    REAL(DP), DIMENSION(2) :: wlim_LH = (/ 19.20, 37.10 /)
  END TYPE spect_ran
  PUBLIC :: spect_ran

  TYPE(spect_ran), SAVE, PUBLIC :: ran

  !!--------------------
  !! Calibration errors
  !!--------------------
  TYPE spect_err
    REAL(DP) :: caliberr_AKARI_NG = 0.1_DP ! Ohyama+2007
    REAL(DP) :: caliberr_SL = 0.05_DP
    REAL(DP) :: caliberr_LL = 0.05_DP
    REAL(DP) :: caliberr_SH = 0.2_DP
    REAL(DP) :: caliberr_LH = 0.2_DP
  END TYPE spect_err
  PUBLIC :: spect_err

  TYPE(spect_err), SAVE, PUBLIC :: err

  !!----------
  !! Continua
  !!----------
  ! 'a-C_man20nm_big      '
  ! 'a-C_man20nm_J13      '
  ! 'a-C_man20nm_small    '
  ! 'a-Enst_Fe_man5nm_J13 '
  ! 'a-Enst_Fe_man5nm     '
  ! 'a-Forst_Fe_man5nm_J13'
  ! 'a-Forst_Fe_man5nm    '
  ! 'ACAR_Z96             ' 
  ! 'ACH2_Z96             ' !!!
  ! 'BE_Z96               ' 
  ! 'FeO_H95              '
  ! 'Gra_D03              '
  ! 'Gra_LD93             '
  ! 'PAHi_DL07_C11        '
  ! 'PAHi_DL07_G11        '
  ! 'PAHi_DL07            '
  ! 'PAHi_LD01            '
  ! 'PAHn_DL07_C11        '
  ! 'PAHn_DL07_G11        '
  ! 'PAHn_DL07            '
  ! 'PAHn_LD01            '
  ! 'SiC_LD93             '
  ! 'Sil_D03              ' !!!
  ! 'Sil_LD01             '
  ! 'Sil_LD93             '
  ! 'Sil_Mg07_J03         '
  ! 'Sil_Mg10_J03         '
  ! 'Sil_Mg15_J03         '
  ! 'Sil_Mg20_J03         '
  ! 'Sil_Mg24_J03         '
  ! 'Sil_WD01             '

  !!-------
  !! Lines
  !!-------
  !! Bernard-Salas et al. (2001)
  TYPE line_type
    CHARACTER(31) :: name = "                              "
    CHARACTER(31) :: label = "                              "
    REAL(DP) :: wave
  END TYPE line_type
  PUBLIC :: line_type
  
  TYPE(line_type), DIMENSION(Nline_max), SAVE, PUBLIC :: TABLine = &
    (/ line_type(name='Bracket alpha', &
                label='Bra', wave=4.052_DP), &
       line_type(name='H!E !NI 6-10', &
                label='HI6-10', wave=5.128657_DP), & ! AKARI_NG (2)
       line_type(name='H!D2!N 0-0 S(7)', &
                label='H2S7', wave=5.51116_DP), &
       line_type(name='[Mg!E !NV] 3P!U2!N-3P!U1!N', &
                label='MgV1', wave=5.6098_DP), &
       line_type(name='Humphreys gamma', &
                label='Huc', wave=5.908213_DP), &
       line_type(name='[K!E !NIV]', &
                label='KIV', wave=5.981_DP), &
       line_type(name='H!D2!N 0-0 S(6)', &
                label='H2S6', wave=6.10856_DP), &
       line_type(name='[Cl!E !NV] 2P!U0!N-2P!U0!N', &
                label='ClV', wave=6.709_DP), &
       line_type(name='H!D2!N 0-0 S(5)', &
                label='H2S5', wave=6.90952_DP), &
       line_type(name='HeII 8-9', &
                label='HeII1', wave=6.947984_DP), &
       line_type(name='[Ar!E !NII] 2P!D3/2!N-2P!D1/2!N', &
                label='ArII', wave=6.985274_DP), &
       line_type(name='[Na!E !NIII] 2P!U0!N-2P!U0!N', &
                label='NaIII', wave=7.3178_DP), &
       line_type(name='Pfund alpha', &
                label='Pfa', wave=7.459858_DP), &
       line_type(name='Humphreys beta', &
                label='Hub', wave=7.502493_DP), & ! SL2 (12)
       line_type(name='[Ne!E !NVI] 2P!U0!N-2P!U0!N', &
                label='NeVI', wave=7.6524_DP), &
       line_type(name='[Fe!E !NVII] 3F!U3!N-3F!U4!N', &
                label='FeVII', wave=7.8145_DP), &
       line_type(name='[Ar!E !NV] 3P!U1!N-3P!U2!N', &
                label='ArV', wave=7.90158_DP), &
       line_type(name='H!D2!N 0-0 S(4)', &
                label='H2S4', wave=8.02505_DP), &
       line_type(name='H!E !NI 7-10', &
                label='HI7-10', wave=8.760064_DP), &
       line_type(name='[Ar!E !NIII] 3P!D2!N-3P!D1!N', &
                label='ArIII1', wave=8.99103_DP), &
       line_type(name='[Ni!E !NVI] 4P!U5/2!N-4F!U5/2!N', &
                label='NiVI', wave=9.042_DP), &
       line_type(name='[Fe!E !NVII] 3F!D2!N-3F!D3!N', &
                label='FeVII', wave=9.5267_DP), &
       line_type(name='H!D2!N 0-0 S(3)', &
                label='H2S3', wave=9.66491_DP), &
       line_type(name='He!E !NII 9-10', &
                label='HeII2', wave=9.713475_DP), &
       line_type(name='[Si!E !NI] 1P!D1!N-1P!D2!N', &
                label='SiI', wave=10.3385_DP), &
       line_type(name='[S!E !NIV] 2P!D3/2!N-2P!D1/2!N', &
                label='SIV', wave=10.5105_DP), &
       line_type(name='H!D2!N 0-0 S(2)', &
                label='H2S2', wave=12.27861_DP), &
       line_type(name='Humphreys alpha', &
                label='Hua', wave=12.368527_DP), &
       line_type(name='[Ne!E !NII] 2P!D3/2!N-2P!D1/2!N', &
                label='NeII', wave=12.81355_DP), &
       line_type(name='[Ar!E !NV] 3P!D0!N-3P!D1!N', &
                label='ArV', wave=13.10219_DP), &
       line_type(name='[Mg!E !NV] 3P!D1!N-3P!D0!N', &
                label='MgV2', wave=13.521_DP), & ! SL1 (17)
       line_type(name='[Ne!E !NV] 3P!D1!N-3P!D2!N', &
                label='NeV1', wave=14.32168_DP), &
       line_type(name='[Ne!E !NIII] 3P!D2!N-3P!D1!N', &
                label='NeIII1', wave=15.555_DP), &
       line_type(name='H!D2!N 0-0 S(1)', &
                label='H2S1', wave=17.03483_DP), &
       line_type(name='H!E !NI 11-18', &
                label='HI11-18', wave=17.608246_DP), &
       line_type(name='[S!E !NIII] 3P!D2!N-3P!D1!N', &
                label='SIII1', wave=18.7129_DP), &
       line_type(name='H!E !NI 7-8', &
                label='HI7-8', wave=19.061898_DP), &
       line_type(name='[Ar!E !NIII] 3P!U1!N-3P!U0!N', &
                label='ArIII2', wave=21.8291_DP), &  ! LL2 (6)
       line_type(name='[Ne!E !NV] 3P!U0!N-3P!U1!N', &
                label='NeV2', wave=24.3175_DP), &
       line_type(name='[O!E !NIV] 2P!D3/2!N-2P!D1/2!N', &
                label='OIV', wave=25.8903_DP), &
       line_type(name='[Fe!E !NII] 6D!D7/2!N-6D!D9/2!N', &
                label='FeII1', wave=25.9882_DP), &
       line_type(name='H!D2!N 0-0 S(0)', &
                label='H2S0', wave=28.21883_DP), &
       line_type(name='[S!E !NIII] 3P!D1!N-3P!D0!N', &
                label='SIII2', wave=33.4810_DP), &
       line_type(name='[Si!E !NII] 2P!D3/2!N-2P!D1/2!N', &
                label='SiII', wave=34.8152_DP), &
       line_type(name='[Fe!E !NII] 6D!D5/2!N-6D!D7/2!N', &
                label='FeII2', wave=35.3491_DP), &
       line_type(name='[Ne!E !NIII] 3P!D1!N-3P!D0!N', &
                label='NeIII2', wave=36.0135_DP) /) ! LL1 (9)

  !!-------
  !! Bands
  !!-------
  TYPE band_type
    CHARACTER(30) :: label = "                              "
    REAL(DP) :: wave, sigmaS, sigmaL
  END TYPE band_type
  PUBLIC :: band_type
  
  TYPE(band_type), DIMENSION(Nband_max), SAVE, PUBLIC :: TABand = &
    (/ band_type(label='Main 3.3', wave=3.291_DP, &
                sigmaS=0.020_DP, sigmaL=0.019_DP), &
       band_type(label='Main 3.4', wave=3.399_DP, &
                sigmaS=0.011_DP, sigmaL=0.024_DP), &
       band_type(label='Small 3.5', wave=3.499_DP, &
                sigmaS=0.077_DP, sigmaL=0.071_DP), & ! AKARI_NG (3)
       band_type(label='Small 5.2', wave=5.2394667_DP, &
                sigmaS=0.025218240_DP, sigmaL=0.058333420_DP), &
       band_type(label='Small 5.7 (1)', wave=5.6437461_DP, &
                sigmaS=0.040000000_DP, sigmaL=0.080000000_DP), &
       band_type(label='Small 5.7 (2)', wave=5.7490305_DP, &
                sigmaS=0.040000000_DP, sigmaL=0.080000000_DP), &
       band_type(label='Small 6.0', wave=6.0106598_DP, &
                sigmaS=0.040000000_DP, sigmaL=0.066646401_DP), &
       band_type(label='Main 6.2 (1)', wave=6.2034273_DP, & ! 8
                sigmaS=0.031300317_DP, sigmaL=0.060000000_DP), &
       band_type(label='Main 6.2 (2)', wave=6.2672596_DP, &
                sigmaS=0.036922155_DP, sigmaL=0.11633640_DP), &
       band_type(label='Small 6.6', wave=6.6273788_DP, &
                sigmaS=0.12000000_DP, sigmaL=0.12000000_DP), &
       band_type(label='Small 6.8', wave=6.8548833_DP, &
                sigmaS=0.080000000_DP, sigmaL=0.080000000_DP), &
       band_type(label='Small 7.1', wave=7.0791725_DP, &
                sigmaS=0.080000000_DP, sigmaL=0.080000000_DP), & ! SL2 (9)
       band_type(label='Plateau 7.7', wave=7.6000000_DP, &
                sigmaS=0.48000000_DP, sigmaL=0.50247515_DP), &
       band_type(label='Main 7.7 (1)', wave=7.6171123_DP, & ! 14
                sigmaS=0.11856752_DP, sigmaL=0.14531480_DP), &
       band_type(label='Main 7.7 (2)', wave=7.8704769_DP, &
                sigmaS=0.16998357_DP, sigmaL=0.24523967_DP), &
       band_type(label='Small 8.3', wave=8.3623706_DP, &
                sigmaS=0.016256724_DP, sigmaL=0.016256724_DP), &
       band_type(label='Main 8.6', wave=8.6204540_DP, & ! 17
                sigmaS=0.18340577_DP, sigmaL=0.13337054_DP), &
       band_type(label='Small 9.5', wave=9.5244838_DP, &
                sigmaS=0.10766965_DP, sigmaL=0.60000000_DP), &
       band_type(label='Small 10.7', wave=10.707132_DP, &
                sigmaS=0.10000000_DP, sigmaL=0.10000000_DP), &
       band_type(label='Small 11.0', wave=11.038349_DP, &
                sigmaS=0.026989462_DP, sigmaL=0.073146141_DP), &
       band_type(label='Main 11.2', wave=11.237893_DP, & ! 21
                sigmaS=0.053507232_DP, sigmaL=0.15254629_DP), &
       band_type(label='Plateau 11.3', wave=11.400432_DP, &
                sigmaS=0.72000000_DP, sigmaL=0.63657944_DP), &
       band_type(label='Small 11.8', wave=11.796389_DP, &
                sigmaS=0.020813349_DP, sigmaL=0.020813349_DP), &
       band_type(label='Small 11.9', wave=11.949674_DP, &
                sigmaS=0.080352283_DP, sigmaL=0.22192473_DP), &
       band_type(label='Main 12.7 (1)', wave=12.626842_DP, &
                sigmaS=0.20000000_DP, sigmaL=0.094424464_DP), &
       band_type(label='Main 12.7 (2)', wave=12.760273_DP, &
                sigmaS=0.080436118_DP, sigmaL=0.14000000_DP), &
       band_type(label='Small 13.6', wave=13.559342_DP, &
                sigmaS=0.15954880_DP, sigmaL=0.16054015_DP), &
       band_type(label='Small 14.2', wave=14.257133_DP, &
                sigmaS=0.15208135_DP, sigmaL=0.058951597_DP), & ! SL1 (16)
       band_type(label='Small 15.6', wave=15.893117_DP, &
                sigmaS=0.17857214_DP, sigmaL=0.20000000_DP), &
       band_type(label='Small 16.4', wave=16.482868_DP, &
                sigmaS=0.10000000_DP, sigmaL=0.058462024_DP), &
       band_type(label='Plateau 17.0', wave=17.082868_DP, &
                sigmaS=0.49775906_DP, sigmaL=0.56119197_DP), &
       band_type(label='Small 17.4', wave=17.428485_DP, &
                sigmaS=0.10000000_DP, sigmaL=0.10000000_DP), &
       band_type(label='Small 17.8', wave=17.771096_DP, &
                sigmaS=0.030799172_DP, sigmaL=0.075249330_DP), &
       band_type(label='Small 18.9', wave=18.925630_DP, &
                sigmaS=0.034553879_DP, sigmaL=0.11570587_DP) /) ! LL2 (6)

END MODULE auxil
