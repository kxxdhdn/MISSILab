MODULE factable

  USE utilities, ONLY: DP
  IMPLICIT NONE
  PRIVATE

  !!-------------------------
  !! Instrumental resolution
  !!-------------------------
  TYPE instr_RES
    !! dlam/lam = dnu/nu
    REAL(DP) :: CAM = 0.010373835_DP
    REAL(DP) :: SL = 0.0055722841_DP
    REAL(DP) :: SH = 0.00075127093_DP
    REAL(DP) :: LL = 0.0057517091_DP
    REAL(DP) :: LH = 0.00070906159_DP
    REAL(DP) :: SWS = 0.0011044189_DP
    REAL(DP) :: SWSfine = 0.00036671469_DP
    REAL(DP) :: AKARI_Ns = 0.00356688_DP

  END TYPE instr_RES
  PUBLIC :: instr_RES
  TYPE(instr_RES), SAVE, PUBLIC :: RES

  !!-----------
  !! Continuum
  !!-----------
  ! 'a-C_man20nm_big      '
  ! 'a-C_man20nm_J13      '
  ! 'a-C_man20nm_small    '
  ! 'a-Enst_Fe_man5nm_J13 '
  ! 'a-Enst_Fe_man5nm     '
  ! 'a-Forst_Fe_man5nm_J13'
  ! 'a-Forst_Fe_man5nm    '
  ! 'ACAR_Z96             '
  ! 'ACH2_Z96             '
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
  ! 'Sil_D03              '
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
  TYPE lines_IN

    CHARACTER(30) :: name = "                              "
    CHARACTER(30) :: label = "                              "
    REAL(DP) :: wave

  END TYPE lines_IN
  PUBLIC :: lines_IN
  TYPE(lines_IN), DIMENSION(46), SAVE, PUBLIC :: LIN = &
    [ lines_IN(name='Bracket alpha', &
        label='Bra', wave=4.052_DP), &
      lines_IN(name='H!E !NI 6-10', &
        label='HI6-10', wave=5.128657_DP), &
      lines_IN(name='H!D2!N 0-0 S(7)', &
        label='H2S7', wave=5.51116_DP), &
      lines_IN(name='[Mg!E !NV] 3P!U2!N-3P!U1!N', &
        label='MgV1', wave=5.6098_DP), &
      lines_IN(name='Humphreys gamma', &
        label='Huc', wave=5.908213_DP), &
      lines_IN(name='[K!E !NIV]', &
        label='KIV', wave=5.981_DP), &
      lines_IN(name='H!D2!N 0-0 S(6)', &
        label='H2S6', wave=6.10856_DP), &
      lines_IN(name='[Cl!E !NV] 2P!U0!N-2P!U0!N', &
        label='ClV', wave=6.709_DP), &
      lines_IN(name='H!D2!N 0-0 S(5)', &
        label='H2S5', wave=6.90952_DP), &
      lines_IN(name='HeII 8-9', &
        label='HeII1', wave=6.947984_DP), &
      lines_IN(name='[Ar!E !NII] 2P!D3/2!N-2P!D1/2!N', &
        label='ArII', wave=6.985274_DP), &
      lines_IN(name='[Na!E !NIII] 2P!U0!N-2P!U0!N', &
        label='NaIII', wave=7.3178_DP), &
      lines_IN(name='Pfund alpha', &
        label='Pfa', wave=7.459858_DP), &
      lines_IN(name='Humphreys beta', &
        label='Hub', wave=7.502493_DP), &
      lines_IN(name='[Ne!E !NVI] 2P!U0!N-2P!U0!N', &
        label='NeVI', wave=7.6524_DP), &
      lines_IN(name='[Fe!E !NVII] 3F!U3!N-3F!U4!N', &
        label='FeVII', wave=7.8145_DP), &
      lines_IN(name='[Ar!E !NV] 3P!U1!N-3P!U2!N', &
        label='ArV', wave=7.90158_DP), &
      lines_IN(name='H!D2!N 0-0 S(4)', &
        label='H2S4', wave=8.02505_DP), &
      lines_IN(name='H!E !NI 7-10', &
        label='HI7-10', wave=8.760064_DP), &
      lines_IN(name='[Ar!E !NIII] 3P!D2!N-3P!D1!N', &
        label='ArIII1', wave=8.99103_DP), &
      lines_IN(name='[Ni!E !NVI] 4P!U5/2!N-4F!U5/2!N', &
        label='NiVI', wave=9.042_DP), &
      lines_IN(name='[Fe!E !NVII] 3F!D2!N-3F!D3!N', &
        label='FeVII', wave=9.5267_DP), &
      lines_IN(name='H!D2!N 0-0 S(3)', &
        label='H2S3', wave=9.66491_DP), &
      lines_IN(name='He!E !NII 9-10', &
        label='HeII2', wave=9.713475_DP), &
      lines_IN(name='[Si!E !NI] 1P!D1!N-1P!D2!N', &
        label='SiI', wave=10.3385_DP), &
      lines_IN(name='[S!E !NIV] 2P!D3/2!N-2P!D1/2!N', &
        label='SIV', wave=10.5105_DP), &
      lines_IN(name='H!D2!N 0-0 S(2)', &
        label='H2S2', wave=12.27861_DP), &
      lines_IN(name='Humphreys alpha', &
        label='Hua', wave=12.368527_DP), &
      lines_IN(name='[Ne!E !NII] 2P!D3/2!N-2P!D1/2!N', &
        label='NeII', wave=12.81355_DP), &
      lines_IN(name='[Ar!E !NV] 3P!D0!N-3P!D1!N', &
        label='ArV', wave=13.10219_DP), &
      lines_IN(name='[Mg!E !NV] 3P!D1!N-3P!D0!N', &
        label='MgV2', wave=13.521_DP), &
      lines_IN(name='[Ne!E !NV] 3P!D1!N-3P!D2!N', &
        label='NeV1', wave=14.32168_DP), &
      lines_IN(name='[Ne!E !NIII] 3P!D2!N-3P!D1!N', &
        label='NeIII1', wave=15.555_DP), &
      lines_IN(name='H!D2!N 0-0 S(1)', &
        label='H2S1', wave=17.03483_DP), &
      lines_IN(name='H!E !NI 11-18', &
        label='HI11-18', wave=17.608246_DP), &
      lines_IN(name='[S!E !NIII] 3P!D2!N-3P!D1!N', &
        label='SIII1', wave=18.7129_DP), &
      lines_IN(name='H!E !NI 7-8', &
        label='HI7-8', wave=19.061898_DP), &
      lines_IN(name='[Ar!E !NIII] 3P!U1!N-3P!U0!N', &
        label='ArIII2', wave=21.8291_DP), &
      lines_IN(name='[Ne!E !NV] 3P!U0!N-3P!U1!N', &
        label='NeV2', wave=24.3175_DP), &
      lines_IN(name='[O!E !NIV] 2P!D3/2!N-2P!D1/2!N', &
        label='OIV', wave=25.8903_DP), &
      lines_IN(name='[Fe!E !NII] 6D!D7/2!N-6D!D9/2!N', &
        label='FeII1', wave=25.9882_DP), &
      lines_IN(name='H!D2!N 0-0 S(0)', &
        label='H2S0', wave=28.21883_DP), &
      lines_IN(name='[S!E !NIII] 3P!D1!N-3P!D0!N', &
        label='SIII2', wave=33.4810_DP), &
      lines_IN(name='[Si!E !NII] 2P!D3/2!N-2P!D1/2!N', &
        label='SiII', wave=34.8152_DP), &
      lines_IN(name='[Fe!E !NII] 6D!D5/2!N-6D!D7/2!N', &
        label='FeII2', wave=35.3491_DP), &
      lines_IN(name='[Ne!E !NIII] 3P!D1!N-3P!D0!N', &
        label='NeIII2', wave=36.0135_DP) ]

  !!-------
  !! Bands
  !!-------
  TYPE bands_IN

    CHARACTER(30) :: label = "                              "
    REAL(DP) :: wave, sigmaS, sigmaL

  END TYPE bands_IN
  PUBLIC :: bands_IN
  TYPE(bands_IN), DIMENSION(33), SAVE, PUBLIC :: BIN = &
    [ bands_IN(label='Main 3.3', wave=3.3_DP, &
        sigmaS=0.04_DP, sigmaL=0.04_DP), &
      bands_IN(label='Main 3.4', wave=3.45, &
        sigmaS=0.04_DP, sigmaL=0.04_DP), &
      bands_IN(label='Small 5.2', wave=5.2394667_DP, &
        sigmaS=0.025218240_DP, sigmaL=0.058333420_DP), &
      bands_IN(label='Small 5.7 (1)', wave=5.6437461_DP, &
        sigmaS=0.040000000_DP, sigmaL=0.080000000_DP), &
      bands_IN(label='Small 5.7 (2)', wave=5.7490305_DP, &
        sigmaS=0.040000000_DP, sigmaL=0.080000000_DP), &
      bands_IN(label='Small 6.0', wave=6.0106598_DP, &
        sigmaS=0.040000000_DP, sigmaL=0.066646401_DP), &
      bands_IN(label='Main 6.2 (1)', wave=6.2034273_DP, &
        sigmaS=0.031300317_DP, sigmaL=0.060000000_DP), &
      bands_IN(label='Main 6.2 (2)', wave=6.2672596_DP, &
        sigmaS=0.036922155_DP, sigmaL=0.11633640_DP), &
      bands_IN(label='Small 6.6', wave=6.6273788_DP, &
        sigmaS=0.12000000_DP, sigmaL=0.12000000_DP), &
      bands_IN(label='Small 6.8', wave=6.8548833_DP, &
        sigmaS=0.080000000_DP, sigmaL=0.080000000_DP), &
      bands_IN(label='Small 7.1', wave=7.0791725_DP, &
        sigmaS=0.080000000_DP, sigmaL=0.080000000_DP), &
      bands_IN(label='Plateau 7.7', wave=7.6000000_DP, &
        sigmaS=0.48000000_DP, sigmaL=0.50247515_DP), &
      bands_IN(label='Main 7.7 (1)', wave=7.6171123_DP, &
        sigmaS=0.11856752_DP, sigmaL=0.14531480_DP), &
      bands_IN(label='Main 7.7 (2)', wave=7.8704769_DP, &
        sigmaS=0.16998357_DP, sigmaL=0.24523967_DP), &
      bands_IN(label='Small 8.3', wave=8.3623706_DP, &
        sigmaS=0.016256724_DP, sigmaL=0.016256724_DP), &
      bands_IN(label='Main 8.6', wave=8.6204540_DP, &
        sigmaS=0.18340577_DP, sigmaL=0.13337054_DP), &
      bands_IN(label='Small 9.5', wave=9.5244838_DP, &
        sigmaS=0.10766965_DP, sigmaL=0.60000000_DP), &
      bands_IN(label='Small 10.7', wave=10.707132_DP, &
        sigmaS=0.10000000_DP, sigmaL=0.10000000_DP), &
      bands_IN(label='Small 11.0', wave=11.038349_DP, &
        sigmaS=0.026989462_DP, sigmaL=0.073146141_DP), &
      bands_IN(label='Main 11.2', wave=11.237893_DP, &
        sigmaS=0.053507232_DP, sigmaL=0.15254629_DP), &
      bands_IN(label='Plateau 11.3', wave=11.400432_DP, &
        sigmaS=0.72000000_DP, sigmaL=0.63657944_DP), &
      bands_IN(label='Small 11.8', wave=11.796389_DP, &
        sigmaS=0.020813349_DP, sigmaL=0.020813349_DP), &
      bands_IN(label='Small 11.9', wave=11.949674_DP, &
        sigmaS=0.080352283_DP, sigmaL=0.22192473_DP), &
      bands_IN(label='Main 12.7 (1)', wave=12.626842_DP, &
        sigmaS=0.20000000_DP, sigmaL=0.094424464_DP), &
      bands_IN(label='Main 12.7 (2)', wave=12.760273_DP, &
        sigmaS=0.080436118_DP, sigmaL=0.14000000_DP), &
      bands_IN(label='Small 13.6', wave=13.559342_DP, &
        sigmaS=0.15954880_DP, sigmaL=0.16054015_DP), &
      bands_IN(label='Small 14.2', wave=14.257133_DP, &
        sigmaS=0.15208135_DP, sigmaL=0.058951597_DP), &
      bands_IN(label='Small 15.6', wave=15.893117_DP, &
        sigmaS=0.17857214_DP, sigmaL=0.20000000_DP), &
      bands_IN(label='Small 16.4', wave=16.482868_DP, &
        sigmaS=0.10000000_DP, sigmaL=0.058462024_DP), &
      bands_IN(label='Plateau 17.0', wave=17.082868_DP, &
        sigmaS=0.49775906_DP, sigmaL=0.56119197_DP), &
      bands_IN(label='Small 17.4', wave=17.428485_DP, &
        sigmaS=0.10000000_DP, sigmaL=0.10000000_DP), &
      bands_IN(label='Small 17.8', wave=17.771096_DP, &
        sigmaS=0.030799172_DP, sigmaL=0.075249330_DP), &
      bands_IN(label='Small 18.9', wave=18.925630_DP, &
        sigmaS=0.034553879_DP, sigmaL=0.11570587_DP) ]

END MODULE factable
