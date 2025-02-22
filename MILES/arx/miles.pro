;;*****************************************************************************
;;*
;;*                                MILES
;;*                  Mid-Infrared Line Extraction Software
;:*
;;*****************************************************************************


  ;;-------------------------------------------------------------------------
  ;;                 Make the Qabs Structure Once and For All
  ;;-------------------------------------------------------------------------

FUNCTION make_Qabs, waveIN, TYPE=type

  ;; Structure initilisation
  Nw = N_ELEMENTS(waveIN)
  nu = !MKS.clight/!MKS.micron/waveIN
  Nbb = N_ELEMENTS(type)
  rad = 0.01D ;; microns
  units = "w, rad [micron]; nu [Hz]; Qova [m-1]; rho [kg/m3]"
  Qabs_str = REPLICATE({Nw:Nw, w:waveIN, rad:rad, type:type, $
                        Qova:DBLARR(Nw), rho:0.D, units:units},Nbb)
  
  ;; Optical properties
  FOR i=0,Nbb-1 DO BEGIN
    CASE type[i] OF
      "carbon_BE": typenorm = "crb"
      "amorphous carbon": typenorm = "crb"
      "silicate": typenorm = "sil"
      "astronomical silicate": typenorm = "sil"
      "graphite": typenorm = "gra"
      "alox_compact": typenorm = "al2o3_compact"
      "compact corundum": typenorm = "al2o3_compact"
      "alox_porous": typenorm = "al2o3_porous"
      "porous corundum": typenorm = "al2o3_porous"
      "alox_compact_cde": typenorm = "al2o3_compact_cde"
      "compact corundum (CDE)": typenorm = "al2o3_compact_cde"
      "alox_koike": typenorm = "al2o3_koike"
      "Koike's corundum": typenorm = "al2o3_koike"
      "alox_koike_cde": typenorm = "al2o3_koike_cde"
      "Koike's corundum (CDE)": typenorm = "al2o3_koike_cde"
      "normal silicates": typenorm = "sil_norm"
      "oxygen deficient silicates": typenorm = "sil_Odef"
      ELSE: typenorm = type[i]
    ENDCASE
    CASE typenorm OF
      "crb": BEGIN
        RESTORE, !IDLFRED+"Dust/Data/MIEwisc_carbon_BE.idl"
        Qabs = REFORM(Qabs_dust[CLOSEST(asize,rad),*])
        Qabs_str[i].Qova = FINTERPOL(Qabs,wave,waveIN,/POW)/rad/!MKS.micron
        Qabs_str[i].rho = 2.24D3
      END
      "sil": BEGIN
        RESTORE, !IDLFRED+"Dust/Data/MIEwisc_silicate.idl"
        Qabs = REFORM(Qabs_sil[CLOSEST(asize,rad),*])
        Qabs_str[i].Qova = FINTERPOL(Qabs,wave,waveIN,/POW)/rad/!MKS.micron
        Qabs_str[i].rho = 3.5D3
      END
      "gra": BEGIN
        RESTORE, !IDLFRED+"Dust/Data/MIEwisc_graphite.idl"
        Qabs = REFORM(Qabs_grf[CLOSEST(asize,rad),*])
        Qabs_str[i].Qova = FINTERPOL(Qabs,wave,waveIN,/POW)/rad/!MKS.micron
        Qabs_str[i].rho = 2.24D3
      END
      "al2o3_compact": BEGIN
        RESTORE, !IDLFRED+"Dust/Data/Qabs_Al2O3_compact.xdr"
        Qabs = REFORM(Qabs_Al2O3[CLOSEST(asize,rad),*])
        Qabs_str[i].Qova = FINTERPOL(Qabs,wave,waveIN,/POW)/rad/!MKS.micron
        Qabs_str[i].rho = 3.5D3
      END
      "al2o3_porous": BEGIN
        RESTORE, !IDLFRED+"Dust/Data/Qabs_Al2O3_porous.xdr"
        Qabs = REFORM(Qabs_Al2O3[CLOSEST(asize,rad),*])
        Qabs_str[i].Qova = FINTERPOL(Qabs,wave,waveIN,/POW)/rad/!MKS.micron
        Qabs_str[i].rho = 3.5D3
      END
      "al2o3_compact_cde": BEGIN
        RESTORE, !IDLFRED+"Dust/Data/Qabs_Al2O3_compact_cde.xdr"
        Qabs = REFORM(Qabs_Al2O3[CLOSEST(asize,rad),*])
        Qabs_str[i].Qova = FINTERPOL(Qabs,wave,waveIN,/POW)/rad/!MKS.micron
        Qabs_str[i].rho = 3.5D3
      END
      "al2o3_koike_cde": BEGIN
        RESTORE, !IDLFRED+"Dust/Data/Qabs_Al2O3_koike_cde.xdr"
        Qabs = REFORM(Qabs_Al2O3[CLOSEST(asize,rad),*])
        Qabs_str[i].Qova = FINTERPOL(Qabs,wave,waveIN,/POW)/rad/!MKS.micron
        Qabs_str[i].rho = 3.5D3
      END
      "al2o3_koike": BEGIN
        RESTORE, !IDLFRED+"Dust/Data/Qabs_Al2O3_koike.xdr"
        Qabs = REFORM(Qabs_Al2O3[CLOSEST(asize,rad),*])
        Qabs_str[i].Qova = FINTERPOL(Qabs,wave,waveIN,/POW)/rad/!MKS.micron
        Qabs_str[i].rho = 3.5D3
      END
      "sil_norm": BEGIN
        RESTORE, !IDLFRED+"Dust/Data/Qabs_silicate_normal.xdr"
        Qabs = REFORM(Qabs_sil[CLOSEST(asize,rad),*])
        Qabs_str[i].Qova = FINTERPOL(Qabs,wave,waveIN,/POW)/rad/!MKS.micron
        Qabs_str[i].rho = 3.5D3
      END
      "sil_Odef": BEGIN
        RESTORE, !IDLFRED+"Dust/Data/Qabs_silicate_oxygen_deficient.xdr"
        Qabs = REFORM(Qabs_sil[CLOSEST(asize,rad),*])
        Qabs_str[i].Qova = FINTERPOL(Qabs,wave,waveIN,/POW)/rad/!MKS.micron
        Qabs_str[i].rho = 3.5D3
      END
    ENDCASE
  ENDFOR

  RETURN, Qabs_str

END


  ;;-------------------------------------------------------------------------
  ;;               Analytical Functions of the Individual Features
  ;;-------------------------------------------------------------------------

;; 1) N BB
;;--------
FUNCTION modifBB, Q_str, TEMPBB=temp

  ;; Fnu [MJy/sr]
  ;; in 1 Msun/dist^2 [Msun/kpc2]
  Fnu = !MKS.Msun/!MKS.kpc^2 * 3.D*!PI/4./Q_str.rho * Q_str.Qova $
      * ABSOLUTEBB(Q_str.w,temp,/WAVE,/BNU,/MKS) / !MKS.Jy

  RETURN, Fnu

END

  ;;-------------------------------------------------------------------------

;; 2) Gauss profile
;;-----------------
FUNCTION gaussLine, waveIN, wave_ref, sigma

  nuIN = !MKS.clight/!MKS.micron / waveIN

  nucen = !MKS.clight/!MKS.micron / wave_ref
  nusig = ( 1./(wave_ref-sigma) - 1./(wave_ref+sigma) ) $
          * !MKS.clight/!MKS.micron / 2.D

  Fnu0_norm = EXP( - (nuIN-nucen)^2 / (2.D*nusig^2) ) / SQRT(2.*!PI*nusig^2)

  RETURN, Fnu0_norm

END

  ;;-------------------------------------------------------------------------

;; 3) Asymmetric lorentz profile
;;------------------------------
FUNCTION lorentzBand, waveIN, wave_ref, sigmaS, sigmaL

  IF (N_PARAMS() EQ 3) THEN sigmaL = sigmaS

  Nw = N_ELEMENTS(waveIN)
  nuIN = !MKS.clight/!MKS.micron / waveIN

  nucen = !MKS.clight/!MKS.micron / wave_ref
  nusigS = ( 1./(wave_ref-sigmaS) - 1./wave_ref ) * !MKS.clight/!MKS.micron
  nusigL = ( 1./wave_ref - 1./(wave_ref+sigmaL) ) * !MKS.clight/!MKS.micron

  whS = WHERE(waveIN LE wave_ref, NS, COMP=whL, NCOMP=NL)
  Fnu0_norm = DBLARR(Nw)

  AS = 2.D/(1.+nusigL/nusigS)
  AL = AS*nusigL/nusigS
  IF (NS GT 0) THEN $
    Fnu0_norm[whS] = AS * nusigS/!PI / ( (nuIN[whS]-nucen)^2 + nusigS^2 )
  IF (NL GT 0) THEN $
    Fnu0_norm[whL] = AL * nusigL/!PI / ( (nuIN[whL]-nucen)^2 + nusigL^2 )

  RETURN, Fnu0_norm

END


  ;;-------------------------------------------------------------------------
  ;;         Automatize the degradation of the spectral resolution
  ;;-------------------------------------------------------------------------

FUNCTION degradeRes, wc, dw, INSTRUMENT=instr

  ;; Line width fixed by the instrumental resolution
  dw_w_CAM = 0.010373835D
  dw_w_SL = 0.0055722841D
  dw_w_SH = 0.00075127093D
  dw_w_LL = 0.0057517091D
  dw_w_LH = 0.00070906159D
  dw_w_SWS = 0.0011044189D
  dw_w_SWSfine = 0.00036671469D
  dw_w_AKARI_Ns = 0.00356688D
  
  ;; Fix the linewidth
  CASE instr OF
    "SH": sigma = dw_w_SH * wc
    "LH": sigma = dw_w_LH * wc
    "SL": sigma = dw_w_SL * wc
    "LL": sigma = dw_w_LL * wc
    "SL-SH": BEGIN
      IF (wc LT 10.) THEN sigma = dw_w_SL * wc
      IF (wc GE 10.) THEN sigma = dw_w_SH * wc
    END
    "SL-LL": BEGIN
      IF (wc LT 12.) THEN sigma = dw_w_SL * wc
      IF (wc GE 12.) THEN sigma = dw_w_LL * wc
    END
    "Ns-SL": BEGIN
      IF (wc LT 5.) THEN sigma = dw_w_AKARI_Ns * wc
      IF (wc GE 5.) THEN sigma = dw_w_SL * wc
    END
    "CAM": sigma = dw_w_CAM * wc
    "SWS": sigma = dw_w_SWS * wc
    "SWSfine": sigma = dw_w_SWSfine * wc
    ELSE: sigma = dw
  ENDCASE

  RETURN, sigma

END


  ;;-------------------------------------------------------------------------
  ;;                    Create the parameter structure
  ;;-------------------------------------------------------------------------

FUNCTION parstruct, parm, Nbb=Nbb, Nline=Nline, Nband=Nband, NAv=NAv,Nstar=Nstar

  par = { massBB: DBLARR(MAXI(Nbb,1)), $
          tempBB: DBLARR(MAXI(Nbb,1)), $
          Iline: DBLARR(MAXI(Nline,1)), $
          Cline: DBLARR(MAXI(Nline,1)), $
          Wline: DBLARR(MAXI(Nline,1)), $
          Iband: DBLARR(MAXI(Nband,1)), $
          Cband: DBLARR(MAXI(Nband,1)), $
          WSband: DBLARR(MAXI(Nband,1)), $
          WLband: DBLARR(MAXI(Nband,1)), $
          Av: DBLARR(MAXI(NAv,1)), $
          massStar: 0.D }

  ;; 1) Continuum
  i0 = 0
  FOR i=0,Nbb-1 DO BEGIN
    par.massBB[i] = parm[i0+2*i]
    par.tempBB[i] = parm[i0+2*i+1]
  ENDFOR
  ;; 2) Lines
  i0 += 2*Nbb
  FOR i=0,Nline-1 DO BEGIN
    par.Iline[i] = parm[i0+3*i]
    par.Cline[i] = parm[i0+3*i+1]
    par.Wline[i] = parm[i0+3*i+2]
  ENDFOR
  ;; 3) Bands
  i0 += 3*Nline
  FOR i=0,Nband-1 DO BEGIN
    par.Iband[i] = parm[i0+4*i]
    par.Cband[i] = parm[i0+4*i+1]
    par.WSband[i] = parm[i0+4*i+2]
    par.WLband[i] = parm[i0+4*i+3]
  ENDFOR
  ;; 4) Av
  i0 += 4*Nband
  FOR i=0,NAv-1 DO par.Av[i] = parm[i0+i]
  ;; 5) Star
  i0 += NAv
  IF (Nstar GT 0) THEN par.massStar = parm[i0]

  RETURN, par

END


  ;;=========================================================================


  ;;-------------------------------------------------------------------------
  ;;                    Independent Fit of the lines
  ;;-------------------------------------------------------------------------

PRO linealone, wave, parm, FnuOUT, _EXTRA=extra

  ;; Line
  Iline = parm[0]
  Cline = parm[1]
  Wline = parm[2]
  Fnu_line = Iline*GAUSSLINE(wave,Cline,Wline)

  ;; Continuum
  a = parm[3]
  b = parm[4]
  c = parm[5]
  x = (wave - extra.center)/extra.center
  Fnu_cont = a*(1+b*x+c*x^2)

  ;; Total
  FnuOUT = Fnu_cont + Fnu_line

  ;; Display
  IF (extra.noplot LT 1) THEN BEGIN
    wh = WHERE(extra.FnuObs GT 0.)
    ranFnu = [MIN(extra.FnuObs[wh])*0.9,1.05*MAX(extra.FnuObs)]
    PLOTAXIS, MI=extra.milli, XRANGE=extra.ranw, YRANGE=ranFnu, $
              YLOG=extra.logplot, TITLE="!6"+extra.title, $
              XTITLE=!PL.wmic, YTITLE=!PL.Fnu+" [MJy/sr]"
    OPLOT, wave, Fnu_cont, COLOR=!RED
    OPLOT, wave, Fnu_cont+Fnu_line, COLOR=!BLUE
    OPLOT, wave, FnuOUT, LINE=2, COLOR=!DARKGREY
    WONDERRPLOT, extra.wavObs, extra.FnuObs, extra.wavObs*0., extra.dFnuObs, $
                 PSYM="circ", SYMSIZE=0.3*!P.SYMSIZE, BARTHICK=0.3*!P.THICK
  ENDIF

END


  ;;-------------------------------------------------------------------------
  ;;                              Total model
  ;;-------------------------------------------------------------------------

PRO specmodel, wave, parm, FnuOUT, FNU_BB=Fnu_cont_tab, PABS=Pabs, $
               FNU_BAND=Fnu_band_tab, FNU_LINE=Fnu_line, $
               FNU_STAR=Fnu_star, _EXTRA=extra

  ;; Parameters
  par = PARSTRUCT(parm,Nbb=extra.Nbb,Nline=extra.Nline,Nband=extra.Nband, $ 
                       NAv=extra.NAv,Nstar=extra.Nstar)
  w = extra.w
  FnuOUT = DBLARR(extra.Nw)


  ;; Screen extinction
  ;;------------------
  IF (extra.NAv GT 0) THEN BEGIN 
    tau_tab = DBLARR(extra.NAv,extra.Nw)
    FOR i=0,extra.NAv-1 DO BEGIN
      Av = par.Av[i]
      tau_tab[i,*] = Av/1.086*extra.extinct[i].tau0
    ENDFOR
    Pabs = EXP(-TOTAL(tau_tab,1))
  ENDIF ELSE Pabs = REPLICATE(1.D,extra.Nw)


  ;; Continuum
  ;;----------
  IF (extra.Nbb GT 0) THEN BEGIN
    Fnu_cont_tab = DBLARR(extra.Nbb,extra.Nw)
    FOR i=0,extra.Nbb-1 DO $
      Fnu_cont_tab[i,*] =par.massBB[i]*MODIFBB(extra.Qabs[i],TEMP=par.tempBB[i])
  ENDIF ELSE Fnu_cont_tab = DBLARR(1,extra.Nw)
  Fnu_cont = TOTAL(Fnu_cont_tab,1)


  ;; Lines
  ;;------
  IF (extra.Nline GT 0) THEN BEGIN

    ;; Line profile
    Fnu_line = DBLARR(extra.Nw)
    IF (extra.Nline2 GT 0) THEN BEGIN
      Fnu_line_tab = DBLARR(extra.Nline,extra.Nw)
      FOR i=0,extra.Nlinefree-1 DO $
        Fnu_line_tab[extra.whlinefree[i],*] = par.Iline[extra.whlinefree[i]] $
          * GAUSSLINE(w,par.Cline[extra.whlinefree[i]], $
                      par.Wline[extra.whlinefree[i]])
      FOR i=0,extra.Nlinefix-1 DO $
        Fnu_line_tab[extra.whlinefix[i],*] = par.Iline[extra.whlinefix[i]] $
          * extra.lines[extra.whlinefix[i]].Fnu0

      Fnu_line += TOTAL(Fnu_line_tab[extra.whline2,*],1)*Pabs
    ENDIF

    ;; Lines which were computed first are not absorbed!
    IF (extra.Nline1 GT 0) THEN Fnu_line += extra.Fnu0_line

  ENDIF ELSE Fnu_line = DBLARR(extra.Nw)


  ;; Bands
  ;;------
  IF (extra.Nband GT 0) THEN BEGIN
    Fnu_band_tab = DBLARR(extra.Nband,extra.Nw)

    ;; Bands which are free to vary
    FOR i=0,extra.Nbandfree-1 DO $
      Fnu_band_tab[extra.whbandfree[i],*] = par.Iband[extra.whbandfree[i]] $
        * LORENTZBAND(w,par.Cband[extra.whbandfree[i]], $
                        par.WSband[extra.whbandfree[i]], $
                        par.WLband[extra.whbandfree[i]])

    ;; Fixed profiles
    FOR i=0,extra.Nbandfix-1 DO $
      Fnu_band_tab[extra.whbandfix[i],*] = par.Iband[extra.whbandfix[i]] $
        * extra.bands[extra.whbandfix[i]].Fnu0

  ENDIF ELSE Fnu_band_tab = DBLARR(1,extra.Nw)
  Fnu_band = TOTAL(Fnu_band_tab,1)


  ;; Stellar continuum
  ;;------------------
  IF (extra.Nstar GT 0) THEN Fnu_star = par.massStar * extra.Fnu_star0 $
                        ELSE Fnu_star = DBLARR(extra.Nw)


  ;; Final fit & display
  ;;--------------------
  Fnu_tot = (Fnu_cont+Fnu_band+Fnu_star)*Pabs + Fnu_line
  IF (extra.hires NE 1) THEN FnuOUT = Fnu_tot $
                        ELSE FnuOUT = FINTERPOL(Fnu_tot,w,extra.wavobs)

  IF (extra.noplot LE 1) THEN BEGIN

    PLOTAXIS, OO=extra.logplot, MI=extra.milli, $
              XRANGE=extra.ranw, YRANGE=extra.ranFnu, $
              TITLE="!6"+extra.title, $
              XTITLE=!PL.wmic, YTITLE=!PL.Fnu+" [MJy/sr]"

      Nw_thresh = 800
      IF (extra.Nw LT Nw_thresh) THEN BEGIN
        OPLOT, extra.wavOBS, extra.FnuOBS
        WONDERRPLOT, extra.wavOBS, extra.FnuOBS, $
                     extra.wavOBS*0., extra.dFnuOBS, $
                     PSYM="circ", SYMSIZE=0.1*!P.SYMSIZE, $
                     BARTHICK=0.1*!P.THICK, COLOR=!DARKGREY
      ENDIF ELSE BEGIN
        CURVEFILL, extra.wavOBS, extra.FnuOBS-extra.dFnuOBS, $
                   extra.wavOBS, extra.FnuOBS+extra.dFnuOBS, $
                   COLOR=!SHADEGREY, XRANGE=MINMAX(extra.wavOBS)
      ENDELSE

      OPLOT, w, Fnu_star*Pabs, COLOR=!YELLOW
      FOR i=0,extra.Nbb-1 DO BEGIN
        CASE extra.cont[i].type OF
          "sil": col = !ROSE
          "sil_norm": col = !ROSE
          "sil_Odef": col = !ROSE
          "al2o3": col = !RED
          "al2o3_compact": col = !RED
          "al2o3_porous": col = !RED
          "al2o3_compact_cde": col = !RED
          "al2o3_koike_cde": col = !RED
          "al2o3_koike": col = !RED
          ELSE: col = !ORANGE
        ENDCASE
        OPLOT, w, REFORM(Fnu_cont_tab[i,*])*Pabs, COLOR=col
      ENDFOR

      FOR i=0,extra.Nline-1 DO $
        OPLOT, w, (Fnu_star+Fnu_cont)*Pabs+Fnu_line, COLOR=!FORREST

      FOR i=0,extra.Nband-1 DO $
        OPLOT, w, (Fnu_star+Fnu_cont+REFORM(Fnu_band_tab[i,*]))*Pabs, $
               COLOR=!BLUE

      IF (extra.Nw LT Nw_thresh) THEN $
        OPLOT, extra.wavOBS, extra.FnuOBS, PSYM=SYMB("circ"), $
               SYMSIZE=0.1*!P.SYMSIZE $
        ELSE OPLOT, extra.wavOBS, extra.FnuOBS, LINE=2
      OPLOT, w, Fnu_tot, COLOR=!LIGHTGREY, THICK=1.3*!P.THICK

    PLOTAXIS, OO=extra.logplot, MI=extra.milli, $
              XRANGE=extra.ranw, YRANGE=extra.ranFnu, $
              TITLE="!6"+extra.title, $
              XTITLE=!PL.wmic, YTITLE=!PL.Fnu+" [MJy/sr]", $
              /NOERASE

 ENDIF


END


  ;;=========================================================================


PRO MILES, wavobsIN, FnuobsIN, dFnuobsIN, WEIGHTS=weightsIN, MASK=maskIN, $
    INSTRUMENT=instrumentIN, $
    ;; Printings and savings
    COPY=filepsIN, PDF=pdfIN, EPS=epsIN, SAVE=filexdrIN, MODELSAVE=modsaveIN, $
    NOPLOT=noplotIN, LOGPLOT=logplotIN, QUIET=quietIN, MILLI=milliIN, $
    TITLE=titleIN, RANGEWAVELENGTH=rangewaveIN, RANGEFLUX=rangefluxIN, $
    HIGHRESOLUTION=hiresIN, COPLINE=coplineIN, $
    ;; Tolerance and iteration
    FTOL=ftolIN, NITER=NiterIN, RESTART=restartIN, $
    ;; Physical components
    BBTYPE=BBtypeIN, LINES=linesIN, BANDS=bandsIN, EXTINCTION=extinctionIN, $
    OLD_STARS=old_starIN, LINE_FIRST=line_firstIN, $
    BB_STRUCTURE=bb_structIN, BAND_STRUCTURE=band_structIN, $
    ;; BB parameters
    TEMPBB=tempBBIN, FIXTEMPBB=fixtempBBIN, TMINBB=TminBBIN, TMAXBB=TmaxBBIN, $
    ;; Lines
    FIXCENTERLINE=fixcenterlineIN, FIXSIGMALINE=fixsigmalineIN, $
    ;; Bands
    FIXCENTERBAND=fixcenterbandIN, FIXSIGMABAND=fixsigmabandIN, $
    UNTIESIGMABAND=untiesigmabandIN, $
    ;; Extinction
    Av=AvIN, FIXAV=fixAvIN, MAXAV=AvmaxIN, TCO2=TCO2IN, $
    ;; Field stars
    MASSSTAR=massStarIN, FIXMASSSTAR=fixmassStarIN, $
    RANGEMASSSTAR=ranmassStarIN, $
    ;; Outputs
    PARSTR=parstrOUT, ERRPARSTR=errparstrOUT, PARINFO=parinfoOUT, $
    FIT=fitOUT, STATUS=statusOUT, CHI2=chi2OUT, NOX=noX

 
  !EXCEPT = 0


  ;;-------------------------------------------------------------------------
  ;;                            Organize the inputs
  ;;-------------------------------------------------------------------------

  ;; Banners
  ;;--------
  ;; Info
  banner = [STRJOIN(REPLICATE("-",80)), $
            STRJOIN(REPLICATE(" ",15))+"MILES: "+ $
            "Mid-Infrared Line Extraction Software.",$
            "It belongs to the SWING (SoftWares for Interpreting Nebulae and " $
            +"Galaxies) package", $
            STRJOIN(REPLICATE(" ",30))+"(written by F. G.).", $
            STRJOIN(REPLICATE(" ",80))]
  
  ;; Manual
  manual = $
  ["MILES, wavOBS[N], FnuOBS[N][Nx,Ny,N], dFnuOBS[N][Nx,Ny,N], INSTRUMENT=[N],"$
   +" WEIGHTS=[N][Nx,Ny,N],", $
   "   MASK=mask[Nx,Ny],", $
   "        (wave in microns is supposed to be sorted, Fnu in MJy/sr)", $
   "   COPY=fileps, /PDF, /EPS, SAVE=filexdr, /MODELSAVE, /NOPLOT, /LOGPLOT,"$
   +" /QUIET, TITLE=''[Nx,Ny]], ", $
   "   /MILLI, RANGEWAVELENGTH=[,], /HIGHRES, FTOL=1.D-5, NITER=250, "$
   +"/RESTART, /COPLINE", $
   "        (Physical components)", $
   "   BBTYPE={1,['crb','gra','sil_norm','al2o3_koike',etc.]}, ", $
   "   LINES={1,['[ArII]','[SIV]',etc.]},", $
   "   BANDS={1,['PAH6.2','PAH7.4',etc.]}, ", $
   "   EXTINCTION={1,['Dudley','CO2']},", $
   "   LINE_FIRST={1,['[ArII]',etc.]},", $
   "        (BB parameters)", $
   "   TEMPBB=[Nbb], FIXTEMPBB={1,[Nbb]}, TMINBB={T,[Nbb]}, TMAXBB={T,[Nbb]},",$
   "        (Field stars)", $
   "   MASSSTAR=, /FIXMASSSTAR, RANGEMASSSTAR=[,], ", $
   "        (Line parameters)", $
   "   FIXCENTERLINE={1,[Nline]}, FIXSIGMALINE={1,[Nline]}, ", $
   "        (band parameters)", $
   "   FIXCENTERBAND={1,[Nband]}, FIXSIGMABAND={1,[Nband]}, /UNTIESIGMABAND, ",$
   "        (Extinction)", $
   "   Av=, /FIXAV, RANGEAV=[,],", $
   "        (Outputs)", $
   "   PARSTR=, ERRPARSTR=, PARINFO=, FIT=, STATUS=, CHI2=," $
   +" BAND_STRUCTURE=bandIN, BB_STRUCTURE=bbIN, /NOX"]
  IF (N_PARAMS() LE 1) THEN BEGIN
    FOR i=0,N_ELEMENTS(manual)-1 DO PRINT, manual[i]
    RETURN
  ENDIF
  t0 = SYSTIME(1)


  ;; Input parameters
  ;;-----------------
  ;; Observed SED (wave is supposed to be sorted and uniqued)
  wavobs = wavobsIN
  nuobs = !MKS.clight/(wavobs*!MKS.micron)
  IF (NOT KEYWORD_SET(rangewaveIN)) THEN ranw = [ MAXI(1.D,0.9*MIN(wavobs)),$
                                                  MINI(40.D,1.1*MAX(wavobs)) ] $
                                    ELSE ranw = rangewaveIN
  IF (N_PARAMS() EQ 2) THEN dFnuobsIN = 0.05*FnuobsIN
  siz = SIZE(FnuOBSIN)
  CASE siz[0] OF
    1: BEGIN
      Nx = 1
      Ny = 1
      Nobs = siz[1]
    END
    2: BEGIN
      Nx = siz[1] 
      Ny = 1
      Nobs = siz[2]
    END
    3: BEGIN
      Nx = siz[1] 
      Ny = siz[2] 
      Nobs = siz[3]
    END
  ENDCASE
  IF (Nobs NE N_ELEMENTS(wavobs)) THEN BEGIN
    PRINT, "!!! MILES: WAVOBS and FNUOBS have incompatible sizes..."
    RETURN
  ENDIF
  Fnuobs = REFORM(FnuobsIN,Nx,Ny,Nobs)
  dFnuobs = REFORM(dFnuobsIN,Nx,Ny,Nobs)
  IF (NOT KEYWORD_SET(maskIN)) THEN BEGIN
    mask = BYTARR(Nx,Ny)
    wh0 = WHERE(TOTAL(Fnuobs,3) EQ 0.D OR FINITE(TOTAL(Fnuobs,3)) EQ 0,N0)
    IF (N0 GT 0) THEN mask[wh0] = 1
  ENDIF ELSE mask = REFORM(maskIN,Nx,Ny)

  ;; Chi square weighting
  IF (NOT KEYWORD_SET(weightsIN)) THEN weights = 1/dFnuobs^2 $
                                  ELSE weights = REFORM(weightsIN,Nx,Ny,Nobs)

  ;; Instrument
  IF (NOT KEYWORD_SET(instrumentIN)) THEN instrument = REPLICATE("",Nobs) $
    ELSE IF (N_ELEMENTS(instrument) EQ Nobs) THEN instrument = instrumentIN $
            ELSE instrument = REPLICATE(instrumentIN[0],Nobs)

  ;; EPS files
  IF KEYWORD_SET(filepsIN) THEN BEGIN
    IF (N_ELEMENTS(filepsIN) EQ 1) THEN BEGIN
      filepsrad = STRREPLACE(filepsIN,[".eps"],[""])
      IF (Ny EQ 1) THEN $
        IF (Nx EQ 1) THEN fileps = filepsrad+".eps" $
                     ELSE fileps = filepsrad+"_"+PRING(INDGEN(Nx)+1)+".eps" $
      ELSE BEGIN
        fileps = STRARR(Nx,Ny)
        FOR x=0,Nx-1 DO $
          fileps[x,*] = filepsrad+"_"+PRING(x+1)+"_"+PRING(INDGEN(Ny)+1)+".eps"
      ENDELSE
    ENDIF ELSE fileps = filepsIN
  ENDIF ELSE fileps = 0


  ;; Fit control parameters
  ;;-----------------------
  ;; Tolerance, iteration, printing
  IF (NOT KEYWORD_SET(milliIN)) THEN milli = 0 ELSE milli = milliIN
  IF (NOT KEYWORD_SET(titleIN)) THEN title = "" ELSE title = titleIN+" "
  IF (NOT KEYWORD_SET(noplotIN)) THEN noplot = 0 ELSE noplot = noplotIN
  IF (NOT KEYWORD_SET(logplotIN)) THEN logplot = 0 ELSE logplot = logplotIN
  IF (NOT KEYWORD_SET(quietIN)) THEN quiet = 0 ELSE quiet = quietIN
  IF (NOT KEYWORD_SET(ftolIN)) THEN ftol = 1.D-10 ELSE ftol = ftolIN
  IF (NOT KEYWORD_SET(NiterIN)) THEN Niter = 250 ELSE Niter = NiterIN
  IF (NOT KEYWORD_SET(restartIN)) THEN restart = 0 ELSE restart = restartIN
  IF (NOT KEYWORD_SET(hiresIN)) THEN hires = 0 ELSE hires = hiresIN

  ;; Banner
  IF (quiet LE 2) THEN BEGIN
    PRINT, banner
  ENDIF


  ;;-------------------------------------------------------------------------
  ;;               Spectral Feature Properties and Templates
  ;;-------------------------------------------------------------------------

  ;; Model wavelength grid
  ;;----------------------
  IF (NOT KEYWORD_SET(hiresIN)) THEN BEGIN
    wmod = wavOBS
    Nwmod = Nobs
  ENDIF ELSE BEGIN
    Nmin = 10 ;; minimum number of wavelengths
    Nw_mic = 100 ;; number of wavelengths per micron
    N_mic = ranw[1]-ranw[0] ;; number of microns
    Nwmod = MAXI(ROUND(Nw_mic*N_mic),Nmin)
    wmod = [RAMP(Nwmod,ranw[0],ranw[1]),wavOBS]
    ;; Reorder and remove duplicates
    wmod = wmod[UNIQ(wmod)]
    wmod = wmod[SORT(wmod)]
    Nwmod = N_ELEMENTS(wmod)
  ENDELSE
  numod = !MKS.clight/!MKS.micron/wmod
  
  ;; Average interval between two consecutive wavelengths
  dwOBS = (MAX(wavOBS)-MIN(wavOBS)) / Nobs


  ;; 1) Grain continuum
  ;;-------------------
  contstr = {type:"crb", T:50.D, Tmin:10.D, Tmax:3000.D, $
             Mfixed:0, Tfixed:0, Fnu0:DBLARR(Nwmod)}
  IF KEYWORD_SET(BBtypeIN) THEN BEGIN

    ;; Dust composition
    BBtype = BBtypeIN  
    IF (PRING(BBtypeIN[0]) EQ '1') THEN $
      BBtype = ["crb","crb","sil_norm","al2o3_koike"]
    Nbb = N_ELEMENTS(BBtype)
    cont = REPLICATE(contstr,Nbb)
    cont.type = BBtype
    Qabs_str = MAKE_QABS(wmod,TYPE=cont.type)

    ;; Initial temperatures
    IF KEYWORD_SET(tempBBIN) THEN BEGIN
      IF (N_ELEMENTS(tempBBIN) NE Nbb) THEN BEGIN
        PRINT, "!!! MILES: wrong TEMPBB keyword..."
        RETURN
      ENDIF ELSE cont.T = tempBBIN
    ENDIF ELSE BEGIN
      IF (PRING(BBtypeIN[0]) EQ '1') THEN BEGIN
        cont.T = [300.D,80.D,150.D,180.D]
      ENDIF ELSE cont.T = 40.D + DINDGEN(Nbb)*30.D
    ENDELSE

    ;; Temperature variation
    IF KEYWORD_SET(fixtempBBIN) THEN BEGIN
      CASE N_ELEMENTS(fixtempBBIN) OF
        1: cont.Tfixed = REPLICATE(fixtempBBIN,Nbb)
        Nbb: cont.Tfixed = fixtempBBIN
        ELSE: BEGIN
          PRINT, "!!! MILES: wrong FIXTEMPBB keyword..."
          RETURN
        END
      ENDCASE
    ENDIF
    
    ;; Minimum temperature
    IF KEYWORD_SET(TminBBIN) THEN BEGIN
      CASE N_ELEMENTS(TminBBIN) OF
        1: cont.Tmin = REPLICATE(TminBBIN,Nbb)
        Nbb: cont.Tmin = TminBBIN
        ELSE: BEGIN
          PRINT, "!!! MILES: wrong TMINBB keyword..."
          RETURN
        END
      ENDCASE      
    ENDIF

    ;; Maximum temperature
    IF KEYWORD_SET(TmaxBBIN) THEN BEGIN
      CASE N_ELEMENTS(TmaxBBIN) OF
        1: cont.Tmax = REPLICATE(TmaxBBIN,Nbb)
        Nbb: cont.Tmax = TmaxBBIN
        ELSE: BEGIN
          PRINT, "!!! MILES: wrong TMAXBB keyword..."
          RETURN
        END
      ENDCASE      
    ENDIF

    ;; Normalized component
    FOR i=0,Nbb-1 DO cont[i].Fnu0 = MODIFBB(Qabs_str[i],TEMPBB=cont[i].T)
    
  ENDIF ELSE BEGIN

    ;; No continuum
    Nbb = 0 
    BBtype = 0
    cont = 0

  ENDELSE

  ;; Special BB structure
  IF KEYWORD_SET(BB_structIN) THEN BEGIN
    cont = BB_structIN
    Nbb = N_ELEMENTS(cont)
    BBtype = cont.type
    Qabs_str = MAKE_QABS(wmod,TYPE=cont.type)
    FOR i=0,Nbb-1 DO cont[i].Fnu0 = MODIFBB(Qabs_str[i],TEMPBB=cont[i].T)
  ENDIF

 
  ;; 2) Lines
  ;;---------
  ;; Bernard-Salas et al. (2001)
  IF KEYWORD_SET(linesIN) THEN BEGIN
  
    lines = REPLICATE( {wave:0.D, label:"", type:"", sigma:0.D, $
                        Ifixed:0, Wfixed:0, Cfixed:0, degpoly:0, $
                        Wlimits:[0.D,0.D], Climits:[0.D,0.D], $
                        Fnu0:DBLARR(Nwmod)}, 51)

    il = 0
    lines[il].wave = 4.052D    &  lines[il].label = "Bra"
      lines[il].type = "Bracket alpha"   &  il += 1
    lines[il].wave = 5.128657D    &  lines[il].label = "HI6-10"
      lines[il].type = "H!E !NI 6-10"   &  il += 1
    lines[il].wave = 5.51116D    &  lines[il].label = "H2S7"
      lines[il].type = "H!D2!N 0-0 S(7)"   &  il += 1
    lines[il].wave = 5.6098D    &  lines[il].label = "MgV"
      lines[il].type = "[Mg!E !NV] 3P!U2!N-3P!U1!N"   &  il += 1
    lines[il].wave = 5.908213D   &  lines[il].label = "Huc"
      lines[il].type = "Humphreys gamma"   &  il += 1
    lines[il].wave = 5.981D   &  lines[il].label = "KIV"
      lines[il].type = "[K!E !NIV]"   &  il += 1
    lines[il].wave = 6.10856D    &  lines[il].label = "H2S6"
      lines[il].type = "H!D2!N 0-0 S(6)"   &  il += 1
    lines[il].wave = 6.709D    &  lines[il].label = "ClV"
      lines[il].type = "[Cl!E !NV] 2P!U0!N-2P!U0!N"   &  il += 1
    lines[il].wave = 6.90952D    &  lines[il].label = "H2S5"
      lines[il].type = "H!D2!N 0-0 S(5)"   &  il += 1
    lines[il].wave = 6.947984D    &  lines[il].label = "HeII1"
      lines[il].type = "HeII 8-9"   &  il += 1
    lines[il].wave = 6.985274D   &  lines[il].label = "ArII"
      lines[il].type = "[Ar!E !NII] 2P!D3/2!N-2P!D1/2!N"   &  il += 1
    lines[il].wave = 7.3178D    &  lines[il].label = "NaIII"
      lines[il].type = "[Na!E !NIII] 2P!U0!N-2P!U0!N"   &  il += 1
    lines[il].wave = 7.459858D   &  lines[il].label = "Pfa"
      lines[il].type = "Pfund alpha"   &  il += 1
    lines[il].wave = 7.502493D   &  lines[il].label = "Hub"
      lines[il].type = "Humphreys beta"  &  il += 1
    lines[il].wave = 7.6524D   &  lines[il].label = "NeVI"
      lines[il].type = "[Ne!E !NVI] 2P!U0!N-2P!U0!N"  &  il += 1
    lines[il].wave = 7.8145D   &  lines[il].label = "FeVII"
      lines[il].type = "[Fe!E !NVII] 3F!U3!N-3F!U4!N"  &  il += 1
    lines[il].wave = 7.90158D   &  lines[il].label = "ArV"
      lines[il].type = "[Ar!E !NV] 3P!U1!N-3P!U2!N"  &  il += 1
    lines[il].wave = 8.02505D    &  lines[il].label = "H2S4"
      lines[il].type = "H!D2!N 0-0 S(4)"  &  il += 1
    lines[il].wave = 8.760064D    &  lines[il].label = "HI7-10"
      lines[il].type = "H!E !NI 7-10"  &  il += 1
    lines[il].wave = 8.99103D    &  lines[il].label = "ArIII1"
      lines[il].type = "[Ar!E !NIII] 3P!D2!N-3P!D1!N"  &  il += 1
    lines[il].wave = 9.042D    &  lines[il].label = "NiVI"
      lines[il].type = "[Ni!E !NVI] 4P!U5/2!N-4F!U5/2!N"  &  il += 1
    lines[il].wave = 9.5267D    &  lines[il].label = "FeVII"
      lines[il].type = "[Fe!E !NVII] 3F!D2!N-3F!D3!N"  &  il += 1
    lines[il].wave = 9.66491D    &  lines[il].label = "H2S3"
      lines[il].type = "H!D2!N 0-0 S(3)"  &  il += 1
    lines[il].wave = 9.713475D    &  lines[il].label = "HeII2"
      lines[il].type = "He!E !NII 9-10"  &  il += 1
    lines[il].wave = 10.3385D    &  lines[il].label = "SiI"
      lines[il].type = "[Si!E !NI] 1P!D1!N-1P!D2!N"  &  il += 1
    lines[il].wave = 10.5105D    &  lines[il].label = "SIV"
      lines[il].type = "[S!E !NIV] 2P!D3/2!N-2P!D1/2!N"  &  il += 1
    lines[il].wave = 12.27861D   &  lines[il].label = "H2S2"
      lines[il].type = "H!D2!N 0-0 S(2)"  &  il += 1
    lines[il].wave = 12.368527D  &  lines[il].label = "Hua"
      lines[il].type = "Humphreys alpha"  &  il += 1
    lines[il].wave = 12.81355D   &  lines[il].label = "NeII"
      lines[il].type = "[Ne!E !NII] 2P!D3/2!N-2P!D1/2!N"  &  il += 1
    lines[il].wave = 13.10219D   &  lines[il].label = "ArV"
      lines[il].type = "[Ar!E !NV] 3P!D0!N-3P!D1!N"  &  il += 1
    lines[il].wave = 13.521D   &  lines[il].label = "MgV"
      lines[il].type = "[Mg!E !NV] 3P!D1!N-3P!D0!N"  &  il += 1
    lines[il].wave = 14.32168D   &  lines[il].label = "NeV1"
      lines[il].type = "[Ne!E !NV] 3P!D1!N-3P!D2!N"  &  il += 1
    lines[il].wave = 15.555D     &  lines[il].label = "NeIII1"
      lines[il].type = "[Ne!E !NIII] 3P!D2!N-3P!D1!N"  &  il += 1
    lines[il].wave = 17.03483D   &  lines[il].label = "H2S1"
      lines[il].type = "H!D2!N 0-0 S(1)"  &  il += 1
    lines[il].wave = 17.608246D   &  lines[il].label = "HI11-18"
      lines[il].type = "H!E !NI 11-18"  &  il += 1
    lines[il].wave = 18.7129D    &  lines[il].label = "SIII1"
      lines[il].type = "[S!E !NIII] 3P!D2!N-3P!D1!N"  &  il += 1
    lines[il].wave = 19.061898D    &  lines[il].label = "HI7-8"
      lines[il].type = "H!E !NI 7-8"  &  il += 1
    lines[il].wave = 21.8291D    &  lines[il].label = "ArIII2"
      lines[il].type = "[Ar!E !NIII] 3P!U1!N-3P!U0!N"  &  il += 1
    lines[il].wave = 24.3175D    &  lines[il].label = "NeV2"
      lines[il].type = "[Ne!E !NV] 3P!U0!N-3P!U1!N"  &  il += 1
    lines[il].wave = 25.8903D    &  lines[il].label = "OIV"
      lines[il].type = "[O!E !NIV] 2P!D3/2!N-2P!D1/2!N"  &  il += 1
    lines[il].wave = 25.9882D    &  lines[il].label = "FeII1"
      lines[il].type = "[Fe!E !NII] 6D!D7/2!N-6D!D9/2!N"
    lines[il].wave = 28.21883D   &  lines[il].label = "H2S0"
      lines[il].type = "H!D2!N 0-0 S(0)"  &  il += 1
    lines[il].wave = 33.4810D    &  lines[il].label = "SIII2"
      lines[il].type = "[S!E !NIII] 3P!D1!N-3P!D0!N"  &  il += 1
    lines[il].wave = 34.8152D    &  lines[il].label = "SiII"
      lines[il].type = "[Si!E !NII] 2P!D3/2!N-2P!D1/2!N"  &  il += 1
    lines[il].wave = 35.3491D    &  lines[il].label = "FeII2"
      lines[il].type = "[Fe!E !NII] 6D!D5/2!N-6D!D7/2!N"  &  il += 1
    lines[il].wave = 36.0135D    &  lines[il].label = "NeIII2"
      lines[il].type = "[Ne!E !NIII] 3P!D1!N-3P!D0!N"  &  il += 1

    ;; Reduce the wavelength range to the observed wavelengths
    lines = lines[WHERE(lines.wave GT MIN(wavOBS) $
                        AND lines.wave LT MAX(wavOBS))]
    Nline = N_ELEMENTS(lines)
    
    ;; Input line list
    IF (PRING(linesIN[0]) NE '1') THEN BEGIN
      Nline = N_ELEMENTS(linesIN)
      whline = INTARR(Nline)
      FOR i=0,Nline-1 DO whline[i] = WHERE(lines.label EQ linesIN[i])
      lines = lines[whline]
    ENDIF ELSE lines = lines[WHERE(lines.wave GT 0.,Nline)]

    ;; Normalised components
    FOR i=0,Nline-1 DO BEGIN
      lines[i].sigma = DEGRADERES(lines[i].wave,dwOBS, $
                         INSTRUMENT=instrument[CLOSEST(wavOBS,lines[i].wave)])
      lines[i].Wlimits = [0.67D,1.5D]*lines[i].sigma
      lines[i].Climits = lines[i].wave + [-0.5,0.5]*dwOBS
      lines[i].Fnu0 = GAUSSLINE(wmod,lines[i].wave,lines[i].sigma)
    ENDFOR

    ;; Constraints
    IF KEYWORD_SET(fixcenterlineIN) THEN BEGIN
      CASE N_ELEMENTS(fixcenterlineIN) OF
        1: lines.Cfixed = REPLICATE(fixcenterlineIN,Nline)
        Nline: lines.Cfixed = fixcenterlineIN
        ELSE: BEGIN
          PRINT, "!!! MILES: wrong FIXCENTERLINE keyword..."
          RETURN
        END
      ENDCASE      
    ENDIF
    IF KEYWORD_SET(fixsigmalineIN) THEN BEGIN
      CASE N_ELEMENTS(fixsigmalineIN) OF
        1: lines.Wfixed = REPLICATE(fixsigmalineIN,Nline)
        Nline: lines.Wfixed = fixsigmalineIN
        ELSE: BEGIN
          PRINT, "!!! MILES: wrong FIXSIGMALINE keyword..."
          RETURN
        END
      ENDCASE
    ENDIF

    ;; First line fit
    IF KEYWORD_SET(line_firstIN) THEN BEGIN
      IF (PRING(line_firstIN[0]) EQ '1') THEN lines.degpoly = 2 ELSE BEGIN
        FOR i=0,Nline-1 DO BEGIN
          wh = WHERE(lines[i].label EQ line_firstIN.label,N)
          IF (N EQ 1) THEN lines[i].degpoly = 2
        ENDFOR
      ENDELSE
    ENDIF

  ENDIF ELSE BEGIN

    Nline = 0
    lines = 0    

  ENDELSE


  ;; 3) Bands
  ;;---------
  IF KEYWORD_SET(bandsIN) THEN BEGIN

    ;; The limits are relative Wlimits*sigma
    bands = REPLICATE( {wave:0.D, label:'', type:'', sigmaS:0.D, sigmaL:0.D, $
                        Ifixed:0, Wfixed:0, Cfixed:0, tiesigma:1, $
                        Wlimits:[0.5D,2.D], Climits:[0.D,0.D], $
                        Fnu0:DBLARR(Nwmod)}, 33)

    bands[0].wave = 3.3  & bands[0].label = 'Main 3.3'
      bands[0].sigmaS = 0.04  & bands[0].sigmaL = 0.04
    bands[1].wave = 3.45  & bands[1].label = 'Main 3.4'
      bands[1].sigmaS = 0.04  & bands[1].sigmaL = 0.04
    bands[2].wave = 5.2394667  & bands[2].label = 'Small 5.2'
      bands[2].sigmaS = 0.025218240  & bands[2].sigmaL = 0.058333420
    bands[3].wave = 5.6437461  & bands[3].label = 'Small 5.7 (1)'
      bands[3].sigmaS = 0.040000000  & bands[3].sigmaL = 0.080000000
    bands[4].wave = 5.7490305  & bands[4].label = 'Small 5.7 (2)'
      bands[4].sigmaS = 0.040000000  & bands[4].sigmaL = 0.080000000
    bands[5].wave = 6.0106598  & bands[5].label = 'Small 6.0'
      bands[5].sigmaS = 0.040000000  & bands[5].sigmaL = 0.066646401
    bands[6].wave = 6.2034273  & bands[6].label = 'Main 6.2 (1)'
      bands[6].sigmaS = 0.031300317  & bands[6].sigmaL = 0.060000000
    bands[7].wave = 6.2672596  & bands[7].label = 'Main 6.2 (2)'
      bands[7].sigmaS = 0.036922155  & bands[7].sigmaL = 0.11633640
    bands[8].wave = 6.6273788  & bands[8].label = 'Small 6.6'
      bands[8].sigmaS = 0.12000000  & bands[8].sigmaL = 0.12000000
    bands[9].wave = 6.8548833  & bands[9].label = 'Small 6.8'
      bands[9].sigmaS = 0.080000000  & bands[9].sigmaL = 0.080000000
    bands[10].wave = 7.0791725  & bands[10].label = 'Small 7.1'
      bands[10].sigmaS = 0.080000000  & bands[10].sigmaL = 0.080000000
    bands[11].wave = 7.6000000  & bands[11].label = 'Plateau 7.7'
      bands[11].sigmaS = 0.48000000  & bands[11].sigmaL = 0.50247515
    bands[12].wave = 7.6171123  & bands[12].label = 'Main 7.7 (1)'
      bands[12].sigmaS = 0.11856752  & bands[12].sigmaL = 0.14531480
    bands[13].wave = 7.8704769  & bands[13].label = 'Main 7.7 (2)'
      bands[13].sigmaS = 0.16998357  & bands[13].sigmaL = 0.24523967
    bands[14].wave = 8.3623706  & bands[14].label = 'Small 8.3'
      bands[14].sigmaS = 0.016256724  & bands[14].sigmaL = 0.016256724
    bands[15].wave = 8.6204540  & bands[15].label = 'Main 8.6'
      bands[15].sigmaS = 0.18340577  & bands[15].sigmaL = 0.13337054
    bands[16].wave = 9.5244838  & bands[16].label = 'Small 9.5'
      bands[16].sigmaS = 0.10766965  & bands[16].sigmaL = 0.60000000
    bands[17].wave = 10.707132  & bands[17].label = 'Small 10.7'
      bands[17].sigmaS = 0.10000000  & bands[17].sigmaL = 0.10000000
    bands[18].wave = 11.038349  & bands[18].label = 'Small 11.0'
      bands[18].sigmaS = 0.026989462  & bands[18].sigmaL = 0.073146141
    bands[19].wave = 11.237893  & bands[19].label = 'Main 11.2'
      bands[19].sigmaS = 0.053507232  & bands[19].sigmaL = 0.15254629
    bands[20].wave = 11.400432  & bands[20].label = 'Plateau 11.3'
      bands[20].sigmaS = 0.72000000  & bands[20].sigmaL = 0.63657944
    bands[21].wave = 11.796389  & bands[21].label = 'Small 11.8'
      bands[21].sigmaS = 0.020813349  & bands[21].sigmaL = 0.020813349
    bands[22].wave = 11.949674  & bands[22].label = 'Small 11.9'
      bands[22].sigmaS = 0.080352283  & bands[22].sigmaL = 0.22192473
    bands[23].wave = 12.626842  & bands[23].label = 'Main 12.7 (1)'
      bands[23].sigmaS = 0.20000000  & bands[23].sigmaL = 0.094424464
    bands[24].wave = 12.760273  & bands[24].label = 'Main 12.7 (2)'
      bands[24].sigmaS = 0.080436118  & bands[24].sigmaL = 0.14000000
    bands[25].wave = 13.559342  & bands[25].label = 'Small 13.6'
      bands[25].sigmaS = 0.15954880  & bands[25].sigmaL = 0.16054015
    bands[26].wave = 14.257133  & bands[26].label = 'Small 14.2'
      bands[26].sigmaS = 0.15208135  & bands[26].sigmaL = 0.058951597
    bands[27].wave = 15.893117  & bands[27].label = 'Small 15.6'
      bands[27].sigmaS = 0.17857214  & bands[27].sigmaL = 0.20000000
    bands[28].wave = 16.482868  & bands[28].label = 'Small 16.4'
      bands[28].sigmaS = 0.10000000  & bands[28].sigmaL = 0.058462024
    bands[29].wave = 17.082868  & bands[29].label = 'Plateau 17.0'
      bands[29].sigmaS = 0.49775906  & bands[29].sigmaL = 0.56119197
    bands[30].wave = 17.428485  & bands[30].label = 'Small 17.4'
      bands[30].sigmaS = 0.10000000  & bands[30].sigmaL = 0.10000000
    bands[31].wave = 17.771096  & bands[31].label = 'Small 17.8'
      bands[31].sigmaS = 0.030799172  & bands[31].sigmaL = 0.075249330
    bands[32].wave = 18.925630  & bands[32].label = 'Small 18.9'
      bands[32].sigmaS = 0.034553879  & bands[32].sigmaL = 0.11570587

    ;; Input band list
    IF (PRING(bandsIN[0]) NE '1') THEN BEGIN
      Nband = N_ELEMENTS(bandsIN)
      whband = INTARR(Nband)
      FOR i=0,Nband-1 DO whband[i] = WHERE(bands.label EQ bandsIN[i])
      bands = bands[whband]
    ENDIF ELSE bands = bands[WHERE(bands.wave GT 0.,Nband)]

    ;; Normalised components
    FOR i=0,Nband-1 DO BEGIN
      bands[i].sigmaS += DEGRADERES(bands[i].wave,dwOBS, $
                           INSTRUMENT=instrument[CLOSEST(wavOBS,bands[i].wave)])
      bands[i].sigmaL += DEGRADERES(bands[i].wave,dwOBS, $
                           INSTRUMENT=instrument[CLOSEST(wavOBS,bands[i].wave)])
      bands[i].Climits = bands[i].wave + [-1.,1.]*dwOBS
      bands[i].Fnu0 = LORENTZBAND(wmod,bands[i].wave,bands[i].sigmaS, $
                                  bands[i].sigmaL)
    ENDFOR

    ;; Constraints
    IF KEYWORD_SET(fixcenterbandIN) THEN BEGIN
      CASE N_ELEMENTS(fixcenterbandIN) OF
        1: bands.Cfixed = REPLICATE(fixcenterbandIN,Nband)
        Nband: bands.Cfixed = fixcenterbandIN
        ELSE: BEGIN
          PRINT, "!!! MILES: wrong FIXCENTERBAND keyword..."
          RETURN
        END
      ENDCASE      
    ENDIF
    IF KEYWORD_SET(fixsigmabandIN) THEN BEGIN
      CASE N_ELEMENTS(fixsigmabandIN) OF
        1: bands.Wfixed = REPLICATE(fixsigmabandIN,Nband)
        Nband: bands.Wfixed = fixsigmabandIN
        ELSE: BEGIN
          PRINT, "!!! MILES: wrong FIXSIGMABAND keyword..."
          RETURN
        END
      ENDCASE
    ENDIF
    IF KEYWORD_SET(untieWbandIN) THEN $
      IF (N_ELEMENTS(untieWbandIN) EQ Nband) THEN $
        bands.tiesigma = 1-untieWbandIN ELSE bands.tiesigma = 1-untieWbandIN[0]

    ;; Reduce the wavelength range to the observed wavelengths
    bands = bands[WHERE(bands.wave GT MIN(wavOBS) $
                        AND bands.wave LT MAX(wavOBS))]
    Nband = N_ELEMENTS(bands)
    
  ENDIF ELSE BEGIN

    Nband = 0
    bands = 0

  ENDELSE

  ;; Special band structure
  IF KEYWORD_SET(band_structIN) THEN BEGIN
    bands = band_structIN
    Nband = N_ELEMENTS(bands)
    FOR i=0,Nband-1 DO $
      bands[i].Fnu0 = LORENTZBAND(wmod,bands[i].wave,bands[i].sigmaS, $
                                  bands[i].sigmaL)
  ENDIF


  ;; 4) Stellar ISRF
  ;;----------------
  ;; [Lsun/Hz/Msun]
  IF KEYWORD_SET(old_starIN) THEN BEGIN
    Fnu_star0 = ISRFpegase(wmod,AGE=5.D3,/INTERP,/WAVE,/LNU,/MKS,/NOVERB) $
              / !MKS.Lsun
    Nstar = 1
    IF KEYWORD_SET(fixmassStarIN) THEN fixmassStar = fixmassStarIN $
                                  ELSE fixmassStar = 0
    IF KEYWORD_SET(ranmassStarIN) THEN BEGIN
      ranmassStar = ranmassStarIN
      limitedmassStar = [1,1]
    ENDIF ELSE BEGIN
      ranmassStar = [0.D,0.D]
      limitedmassStar = [1,0]
    ENDELSE
  ENDIF ELSE BEGIN
    Nstar = 0
    Fnu_star0 = 0
  ENDELSE


  ;; 5) Dust and ice extinction
  ;;---------------------------
  ;; Constraints
  IF KEYWORD_SET(extinctionIN) THEN BEGIN

    NAv = N_ELEMENTS(extinctionIN)

    IF (PRING(extinctionIN[0]) NE '1') THEN BEGIN
      IF KEYWORD_SET(AvIN) THEN Av = AvIN ELSE Av = REPLICATE(1.D,NAv)
      IF KEYWORD_SET(fixAvIN) THEN $
        IF (N_ELEMENTS(fixAvIN) EQ 1) THEN fixAv = REPLICATE(fixAvIN,NAv) $
                                      ELSE fixAv = fixAvIN $
      ELSE fixAv = REPLICATE(0,NAv)
      IF KEYWORD_SET(AvmaxIN) THEN Avmax = AvmaxIN $
                              ELSE Avmax = REPLICATE(100.D,NAv)
      extinction = extinctionIN
    ENDIF ELSE BEGIN
      IF KEYWORD_SET(AvIN) THEN Av = AvIN ELSE Av = 1.D
      IF KEYWORD_SET(fixAvIN) THEN fixAv = fixAvIN ELSE fixAv = 0
      IF KEYWORD_SET(AvmaxIN) THEN Avmax = AvmaxIN ELSE Avmax = 100.D
      extinction = "dust"
    ENDELSE

  ENDIF ELSE BEGIN

    Av = 0.D   &  fixAv = 1         &   Avmax = 0.D 
    NAv = 0    &  extinction = 0

  ENDELSE

  ;; - a) Dust absorption 
  IF (PRING(extinction[0]) NE '0') THEN BEGIN
    extinction[WHERE(extinction EQ "dust")] = "Dudley"
    ext_law = ["Dudley","GC","Mathis","chiar05","chiar05_GC"]
    FOR i=0,NAv-1 DO BEGIN
      IF (WHERE(ext_law EQ extinction[i]) GE 0) THEN BEGIN
        extfile = !IDLFRED+"Dust/Extinction/law_Av_"+extinction[i]+".dat"
        READCOL, extfile, wave_tab, Av_tab
        ;; Normalization to Av=1
        Av_tab = Av_tab / FINTERPOL(Av_tab,wave_tab,!BAND.V,/POW)
        tau_tab = Av_tab / 1.086
        tau_dust = FINTERPOL(tau_tab,wave_tab,wmod)
        wh0 = WHERE(tau_dust LT 0.,N0)
        IF (N0 GT 0) THEN tau_dust[wh0] = 0.
      ENDIF
    ENDFOR
  ENDIF ELSE tau_dust = DBLARR(Nwmod)

  ;; - b) CO2 ice absorption
  IF KEYWORD_SET(TCO2IN) THEN TCO2 = TCO2IN ELSE TCO2 = 75.D ;; K
  RESTORE, !IDLFRED+"Dust/Data/tau_CO2_all.xdr"
  tau_ice = $
    REFORM(tau_CO2[WHERE(fH2O EQ 1. AND fCH3OH EQ 0.),CLOSEST(T_CO2,TCO2),*])
  whIN = WHERE(ISINE(wmod,MINMAX(wave_CO2)) EQ 1,NIN)
  tau_CO2 = DBLARR(Nwmod)
  IF (NIN GT 0) THEN tau_CO2[whIN] = FINTERPOL(tau_ice,wave_CO2,wmod[whIN])

  ;; All the components
  extinctstr = {w_ref: 0.D, type: "", tau0: DBLARR(Nwmod), fixed: 0, $
                temperature: 0.D, tempfixed: 1, Avmax:0.D }
  IF (NAv GT 0) THEN extinct = REPLICATE(extinctstr,NAv) ELSE extinct = 0
  FOR i=0,NAv-1 DO BEGIN
    IF (extinction[i] EQ "CO2") THEN BEGIN
      extinct[i].w_ref = 15.2D
      extinct[i].tau0 = tau_CO2
      extinct[i].temperature = TCO2
      extinct[i].tempfixed = 1
    ENDIF ELSE BEGIN
      extinct[i].w_ref = !BAND.V 
      extinct[i].tau0 = tau_dust
    ENDELSE
    extinct[i].fixed = fixAV[i]
    extinct[i].type = extinction[i]
    extinct[i].Avmax = Avmax[i]
  ENDFOR


  ;;-------------------------------------------------------------------------
  ;;                Initital Guesses and Parameter Variations
  ;;-------------------------------------------------------------------------

  ;; Parameters
  ;;-----------
  Nparm = 2*Nbb + 3*Nline + 4*Nband + NAv + Nstar
  parinfo = REPLICATE({fixed:0,limited:[0,0],limits:[0.D,0.D],parname:"", $
                       value:0.D,tied:"",mpprint:1},Nx,Ny,Nparm)


  ;; Output arrays
  parstrOUT = REPLICATE(PARSTRUCT(DBLARR(Nparm),Nbb=Nbb,Nline=Nline, $
                        Nband=Nband,NAv=NAv,Nstar=Nstar),Nx,Ny) 
  errparstrOUT = parstrOUT
  statusOUT = INTARR(Nx,Ny)
  chi2OUT = DBLARR(Nx,Ny)
  timin = DBLARR(Nx,Ny)
  IF KEYWORD_SET(modsaveIN) THEN BEGIN
    Fnu_tot = DBLARR(Nx,Ny,Nwmod)
    IF (Nbb GT 0) THEN Fnu_BB = DBLARR(Nx,Ny,Nbb,Nwmod) $
                  ELSE Fnu_BB = DBLARR(Nx,Ny,1,Nwmod)
    Fnu_Line = DBLARR(Nx,Ny,Nwmod)
    IF (Nband GT 0) THEN Fnu_Band = DBLARR(Nx,Ny,Nband,Nwmod) $
                    ELSE Fnu_Band = DBLARR(Nx,Ny,1,Nwmod)
    Pabs = DBLARR(Nx,Ny,Nwmod)
    Fnu_star = DBLARR(Nx,Ny,Nwmod)
  ENDIF
  fitOUT = 0


  ;; BIG LOOP
  ;;---------
FOR x=0,Nx-1 DO BEGIN
FOR y=0,Ny-1 DO BEGIN
IF (mask[x,y] EQ 0) THEN BEGIN
  t1 = SYSTIME(1)

  Fnuobs0 = REFORM(Fnuobs[x,y,*])
  dFnuobs0 = REFORM(dFnuobs[x,y,*])
  weights0 = REFORM(weights[x,y,*])

  IF (quiet LE 2) THEN IF (Nx GT 1) THEN IF (Ny GT 1) THEN $
    PRINT, " - x="+PRING(x+1)+"/"+PRING(Nx)+", y="+PRING(y+1)+"/"+PRING(Ny) $
                                    ELSE PRINT, " - x="+PRING(x+1)+"/"+PRING(Nx)

  ;; Plotting range
  IF KEYWORD_SET(rangefluxIN) THEN ranFnu = rangefluxIN ELSE BEGIN
    wh = WHERE(FnuOBS0-dFnuOBS0 GT 0.,N0)
    IF (N0 GE 10) THEN minOBS = MIN((SMOOTH(FnuOBS0-dFnuOBS0,10))[wh]) $
                  ELSE minOBS = 0.
    maxOBS = ABS(MAX(FnuOBS0+dFnuOBS0))
    IF (minOBS LT 1.D-3*maxOBS) THEN minOBS = 1.D-3*maxOBS
    IF (logplot EQ 0) THEN ranFnu = [0.,1.1*maxOBS] $
                      ELSE ranFnu = [0.5*minOBS,1.1*maxOBS]
  ENDELSE


  ;; Extra parameters
  title_plot = title
  IF (Nx GT 1) THEN $
    IF (Ny GT 1) THEN $
      title_plot += "[x="+PRING(x+1)+"/"+PRING(Nx)+",y=" $
                    +PRING(y+1)+"/"+PRING(Ny)+"]" $
    ELSE title_plot += " [x="+PRING(x+1)+"/"+PRING(Nx)+"]"


  ;;-------------------------------------------------------------------------
  ;;                          Preliminary Line Fit
  ;;-------------------------------------------------------------------------

  IF KEYWORD_SET(line_firstIN) THEN BEGIN
    PRINT, " Preliminary line fit..."
    Iline_first = DBLARR(Nline)
    dIline_first = DBLARR(Nline)
    Cline_first = DBLARR(Nline)
    Wline_first = DBLARR(Nline)
    FOR i=0,Nline-1 DO BEGIN
      IF (lines[i].degpoly GT 0) THEN BEGIN

        ;; Preparation
        ;;------------
        ;; Structure passed by the _EXTRA mechanism
        ranw0 = lines[i].wave + [-6.D,6.D]*lines[i].sigma
        whObs = WHERE(ISINE(wavObs,ranw0),Nobs0)
        whMod = WHERE(ISINE(wmod,ranw0),Nwmod0)
        extra0 = { $
                  ;; Observations
                  Nobs:Nobs0, wavobs:wavobs[whObs], nuobs:nuobs[whObs], $
                  Fnuobs:Fnuobs0[whObs], dFnuobs:dFnuobs0[whObs], $
                  ;; Plot settings
                  title:title_plot+" "+lines[i].label, noplot:noplot, $
                  milli:milli, hires:hires, logplot:logplot, $
                  ;; Fit settings
                  ranw:ranw0, $
                  ;; Fit templates
                  Nw:Nwmod0, w:wmod[whMod], nu:numod[whMod], $
                  center:lines[i].wave } 

        ;; Initial guesses
        Nparm0 = 6
        parm0 = DBLARR(Nparm0)
        parinf0 = REPLICATE(parinfo[0,0,0],Nparm0)
        ;; 1) Lines
        whc = CLOSEST(wavOBS,extra0.center)
        intens = ABS(Fnuobs0[whc])/lines[i].Fnu0[whc] / 2.
        parm0[0:2] = [ intens, lines[i].wave, lines[i].sigma ]
        parinf0[0:2].value = parm0[0:2]
        parinf0[0:2].fixed = [lines[i].Ifixed,lines[i].Cfixed,lines[i].Wfixed]
        parinf0[0].limited = [1,0] 
        parinf0[0].limits = [0.D,0.D]
        parinf0[1].limited = [1,1]
        parinf0[1].limits = lines[i].Climits
        parinf0[2].limited = [1,1]
        parinf0[2].limits = lines[i].Wlimits
        parinf0[0:2].parname = [ "I("+lines[i].label+")", $
                                 "C("+lines[i].label+")", $
                                 "W("+lines[i].label+")" ]
        ;; 2) Grain continuum
        mask0 = (INDGEN(3) LE lines[i].degpoly)
        parm0[3:5] = [ Fnuobs0[whc], 10.D, 10.D ]*mask0
        parinf0[3:5].value = parm0[3:5]
        parinf0[3:5].fixed = 1-mask0
        parinf0[3].limited = [1,0]
        parinf0[3].limits = [0.D,0.D]
        parinf0[3:5].parname = ["a","b","c"]

        ;; Make sure the parameters don't excess the limits
        FOR j=0,Nparm0-1 DO BEGIN
          IF (parinf0[j].fixed EQ 0) THEN BEGIN
            IF (parinf0[j].limited[0] EQ 1 $
                AND parm0[j] LT parinf0[j].limits[0]) THEN BEGIN
              IF (quiet EQ 0) THEN $
                PRINT, "!!! MILES: parameter "+parinf0[j].parname $
                     + " is lower than limit..."
              parm0[j] = parinf0[j].limits[0]
            ENDIF
            IF (parinf0[j].limited[1] EQ 1 $
                AND parm0[j] GT parinf0[j].limits[1]) THEN BEGIN
              IF (quiet EQ 0) THEN $
                PRINT, "!!! MILES: parameter "+parinf0[j].parname $
                     + " is higher than limit..."
              parm0[j] = parinf0[j].limits[1]
            ENDIF
          ENDIF
        ENDFOR

        ;; Printing
        whnoprint = WHERE(parinf0.tied NE "" OR parinf0.fixed EQ 1,Nnoprint)
        IF (Nnoprint GT 0) THEN parinf0[whnoprint].mpprint = 0


        ;; Fit of the lines
        ;;-----------------
        IF (N_ELEMENTS(WHERE(parinf0.fixed EQ 0 $
                             AND parinf0.tied EQ "")) GT Nobs0) THEN BEGIN
          PRINT, "!!! MILES: Number of parameters is higher than number of " $
               + "observations for the separate line fit..."
          RETURN
        ENDIF
        IF (quiet GE 1) THEN quietFIT = 1 ELSE quietFIT = 0

        ;; Actual fit
        weiline = weights0[whObs]
        fit = MPCURVEFIT(extra0.wavObs,extra0.Fnuobs,weiline,parm0, $
                         FUNCTION_NAME="linealone",PARINFO=parinf0, $
                         /NODERIVATIVE,FUNCTARGS=extra0,FTOL=ftol, $
                         XTOL=1.D-50, GTOL=1.D-50,QUIET=quietFIT,COVAR=covar0)

        ;; Error
        sigpar0 = DBLARR(Nparm0)
        FOR j=0,Nparm0-1 DO sigpar0[j] = SQRT(covar0[j,j])

        ;; Final parameters
        Iline_first[i] = parm0[0]
        Cline_first[i] = parm0[1]
        Wline_first[i] = parm0[2]
        dIline_first[i] = sigpar0[0]
        IF (quiet LE 2) THEN $
          PRINT, STRING(lines[i].label,F="(A15)")+"   I="+PRING(Iline_first[i])

        ;; Update the template for the complete fit
        lines[i].Fnu0 = GAUSSLINE(wmod,Cline_first[i],Wline_first[i])

        ;; Save the plot
        IF KEYWORD_SET(coplineIN) THEN BEGIN
          extra0.noplot = 0
          fileps0 = STRREPLACE(fileps[x,y],".eps","_"+lines[i].label+".eps")
          SETPLOT_PS, /TeX, FILE=fileps0, /DEMI
          LINEALONE, extra0.wavObs, parm0, Fnu0, _EXTRA=extra0
          ENDPS, FILE=fileps0, /PDF, NOX=KEYWORD_SET(noX)
        ENDIF

      ENDIF
    ENDFOR
    
    ;; Final template for the line fit
    whtemplate = WHERE(lines.degpoly GT 0,Ntemplate)
    Fnu0_line = DBLARR(Nwmod)
    FOR i=0,Ntemplate-1 DO $
      Fnu0_line += Iline_first[whtemplate[i]]*lines[whtemplate[i]].Fnu0

  ENDIF ELSE Fnu0_line = DBLARR(Nwmod)


  ;;-------------------------------------------------------------------------
  ;;                             Complete Fit
  ;;-------------------------------------------------------------------------

  ;; Initialisation
  ;;---------------
  IF KEYWORD_SET(line_firstIN) THEN PRINT, " Complete fit..."
  parm = DBLARR(Nparm)
  ;; 1) Grain continuum
  i0 = 0
  FOR i=0,Nbb-1 DO BEGIN
    whmax = WHERE(cont[i].Fnu0 EQ MAX(cont[i].Fnu0),Nmax)
    mass = ABS(Fnuobs0[whmax])/cont[i].Fnu0[whmax] / 2.
    i1 = 2*i
    i2 = 2*i+1
    parm[i1:i2] = [mass,cont[i].T]
    parinfo[x,y,i1:i2].value = REFORM(parm[i1:i2],1,1,2)
    parinfo[x,y,i1:i2].fixed = REFORM([cont[i].Mfixed,cont[i].Tfixed],1,1,2)
    parinfo[x,y,i1].limited = [1,0] 
    parinfo[x,y,i1].limits = [0.D,0.D]
    parinfo[x,y,i2].limited = [1,1]
    parinfo[x,y,i2].limits = [cont[i].Tmin,cont[i].Tmax]
    parinfo[x,y,i1:i2].parname = $
      REFORM([ "Mass("+BBtype[i]+"["+PRING(i+1)+"])", $
               "Temp("+BBtype[i]+"["+PRING(i+1)+"])" ],1,1,2)
  ENDFOR
  ;; 2) Lines
  i0 += 2*Nbb
  FOR i=0,Nline-1 DO BEGIN
    i1 = i0+3*i
    i2 = i0+3*i+2
    IF (lines[i].degpoly EQ 0) THEN BEGIN
      whc = CLOSEST(wavOBS,lines[i].wave)
      intens = ABS(Fnuobs0[whc])/lines[i].Fnu0[whc] / 2.
      parm[i1:i2] = [ intens, lines[i].wave, lines[i].sigma ]
      parinfo[x,y,i1:i2].fixed = REFORM([lines[i].Ifixed,lines[i].Cfixed, $
                                         lines[i].Wfixed],1,1,3)
      parinfo[x,y,i1].limited = [1,0] 
      parinfo[x,y,i1].limits = [0.D,0.D]
      parinfo[x,y,i1+1].limited = [1,1]
      parinfo[x,y,i1+1].limits = lines[i].Climits
      parinfo[x,y,i1+2].limited = [1,1]
      parinfo[x,y,i1+2].limits = lines[i].Wlimits
    ENDIF ELSE BEGIN
      whc = CLOSEST(wavOBS,lines[i].wave)
      intens = ABS(Fnuobs0[whc])/lines[i].Fnu0[whc] / 2.
      parm[i1:i2] = [ Iline_first[i], Cline_first[i], Wline_first[i] ]
      parinfo[x,y,i1:i2].fixed = REFORM([1,1,1],1,1,3)
    ENDELSE
    parinfo[x,y,i1:i2].value = REFORM(parm[i1:i2],1,1,3)
    parinfo[x,y,i1:i2].parname = REFORM([ "I("+lines[i].label+")", $
                                          "C("+lines[i].label+")", $
                                          "W("+lines[i].label+")" ],1,1,3)
  ENDFOR
  ;; 3) Bands
  i0 += 3*Nline
  FOR i=0,Nband-1 DO BEGIN
    i1 = i0+4*i
    i2 = i0+4*i+3
    whc = CLOSEST(wavOBS,bands[i].wave)
    intens = ABS(Fnuobs0[whc])/bands[i].Fnu0[whc] / 2.
    parm[i1:i2] = REFORM([ intens, bands[i].wave, bands[i].sigmaS, $
                           bands[i].sigmaL ],1,1,4)
    parinfo[x,y,i1:i2].value = REFORM(parm[i1:i2],1,1,4)
    parinfo[x,y,i1:i2].fixed = REFORM([bands[i].Ifixed,bands[i].Cfixed, $
                                       bands[i].Wfixed,bands[i].Wfixed],1,1,4)
    parinfo[x,y,i1].limited = [1,0] 
    parinfo[x,y,i1].limits = [0.D,0.D]
    parinfo[x,y,i1+1].limited = [1,1]
    parinfo[x,y,i1+1].limits = bands[i].Climits
    parinfo[x,y,i1+2].limited = [1,1]
    parinfo[x,y,i1+2].limits = bands[i].Wlimits*bands[i].sigmaS
    parinfo[x,y,i1+3].limited = [1,1]
    parinfo[x,y,i1+3].limits = bands[i].Wlimits*bands[i].sigmaL
    IF (bands[i].tiesigma) THEN BEGIN
      parinfo[x,y,i1+3].tied = PRING(parm[i1+3]/parm[i1+2]) $
                               +"*P["+PRING(i1+2)+"]"
    ENDIF
    parinfo[x,y,i1:i2].parname = REFORM([ "I("+bands[i].label+")", $
                                          "C("+bands[i].label+")", $
                                          "WS("+bands[i].label+")", $
                                          "WL("+bands[i].label+")" ],1,1,4)
  ENDFOR
  ;; 4) Dust extinction
  i0 += 4*Nband
  FOR i=0,NAv-1 DO BEGIN
    i1 = i0+i
    parm[i1] = 1.D
    parinfo[x,y,i1].value = parm[i1]
    parinfo[x,y,i1].fixed = extinct[i].fixed
    parinfo[x,y,i1].limits = [0.D,extinct[i].Avmax]
    parinfo[x,y,i1].limited = [1,1]
    IF (extinct[i].type EQ "CO2") THEN parinfo[x,y,i1].parname = "A(CO2)" $
                                  ELSE parinfo[x,y,i1].parname = "Av(dust)"
  ENDFOR
  ;; 5) Stellar continuum
  i0 += NAv
  IF (Nstar EQ 1) THEN BEGIN
    whmax = WHERE(Fnu_star0 EQ MAX(Fnu_star0),Nmax)
    mass = ABS(Fnuobs0[whmax])/Fnu_star0[whmax] / 2.
    parm[i0] = mass
    parinfo[x,y,i0].value = mass
    parinfo[x,y,i0].fixed = fixmassStar
    parinfo[x,y,i0].limited = limitedmassStar
    parinfo[x,y,i0].limits = ranmassStar
    parinfo[x,y,i0].parname = "Mass(star)"
  ENDIF

  ;; Make sure the parameters don't excess the limits
  FOR i=0,Nparm-1 DO BEGIN
    IF (parinfo[x,y,i].fixed EQ 0) THEN BEGIN
      IF (parinfo[x,y,i].limited[0] EQ 1 $
          AND parm[i] LT parinfo[x,y,i].limits[0]) THEN BEGIN
        IF (quiet EQ 0) THEN $
          PRINT, "!!! MILES: parameter "+parname[i]+" is lower than limit..."
        parm[i] = parinfo[x,y,i].limits[0]
      ENDIF
      IF (parinfo[x,y,i].limited[1] EQ 1 $
          AND parm[i] GT parinfo[x,y,i].limits[1]) THEN BEGIN
        IF (quiet EQ 0) THEN $
          PRINT, "!!! MILES: parameter "+parname[i]+" is higher than limit..."
        parm[i] = parinfo[x,y,i].limits[1]
      ENDIF
    ENDIF
  ENDFOR

  ;; Printing
  parinfo[x,y,WHERE(parinfo[x,y,*].tied NE "" $
                    OR parinfo[x,y,*].fixed EQ 1)].mpprint = 0

  ;; Separate the treatment of the various lines and bands
  IF (Nline GT 0) THEN BEGIN
    ;; - Lines that were computed previously (for the absorption):
    whline1 = WHERE(lines.degpoly GT 0,Nline1,COMP=whline2,NCOMP=Nline2)
    ;; - Lines where the center and the width are fixed and not
    ;;   previously computed:
    whlinefix = WHERE(lines.degpoly EQ 0 AND $
                      lines.Cfixed EQ 1 AND lines.Wfixed EQ 1 $
                      AND lines.degpoly EQ 0,Nlinefix)
    whlinefree = WHERE(lines.degpoly EQ 0 AND $
                       (lines.Cfixed EQ 0 OR lines.Wfixed EQ 0) $
                       AND lines.degpoly EQ 0,Nlinefree)
  ENDIF ELSE BEGIN
    whline1 = 0
    whline2 = 0 
    whlinefix = 0
    whlinefree = 0
    Nline1 = 0
    Nline2 = 0
    Nlinefix = 0
    Nlinefree = 0
  ENDELSE
  ;; - Bands fixed:
  whbandfix = WHERE(bands.Cfixed EQ 1 AND bands.Wfixed EQ 1,Nbandfix, $
                    COMP=whbandfree,NCOMP=Nbandfree)

  ;; Structure passed by the _EXTRA mechanism
  extra = { $
            ;; Observations
            Nobs:Nobs, wavobs:wavobs, nuobs:nuobs, $
            Fnuobs:Fnuobs0, dFnuobs:dFnuobs0, $
            ;; Plot settings
            title:title_plot, noplot:noplot, milli:milli, hires:hires, $
            ranw:ranw, ranFnu:ranFnu, logplot:logplot, $
            ;; Fit settings
            Nbb:Nbb, Nline:Nline, Nband:Nband, NAv:NAv, Nstar:Nstar, $
            ;; Fit templates
            Nw:Nwmod, w:wmod, nu:numod, cont:cont, lines:lines, bands:bands, $
            extinct:extinct, Fnu_star0:Fnu_star0, Qabs:Qabs_str, $
            Nline1:Nline1, whline1:whline1, Nline2:Nline2, whline2:whline2, $
            Nlinefix:Nlinefix, whlinefix:whlinefix, Nlinefree:Nlinefree, $
            whlinefree:whlinefree, Fnu0_line:Fnu0_line, $
            whbandfix:whbandfix, Nbandfix:Nbandfix, whbandfree:whbandfree, $
            Nbandfree:Nbandfree } 


  ;; Fit of the SED
  ;;---------------
  IF (N_ELEMENTS(WHERE(parinfo[x,y,*].fixed EQ 0 $
                       AND parinfo[x,y,*].tied EQ "")) GT Nobs) THEN BEGIN
    PRINT, "!!! MILES: Number of parameters is higher than number of " $
         + "observations..."
    RETURN
  ENDIF
  IF (quiet GT 1) THEN quietFIT = 1 ELSE quietFIT = 0
  IF (Niter GT 1) THEN BEGIN
    ;; Actual fit
    fit = MPCURVEFIT( wavobs,Fnuobs0,weights0,parm,$
                  FUNCTION_NAME="specmodel",PARINFO=REFORM(parinfo[x,y,*]), $
                  STATUS=status1D,/NODERIVATIVE,FUNCTARGS=extra,FTOL=ftol, $
                  XTOL=1.D-50, GTOL=1.D-50,QUIET=quietFIT,COVAR=covar, $
                  CHI=chi21D,ITMAX=Niter )

    ;; Iteration 
    IF (restart GE 1) THEN BEGIN
pc
stop
      FOR k=0,restart-1 DO BEGIN
        IF (quiet LE 1) THEN PRINT, " Iteration..."
        parinf = REFORM(parinfo[x,y,*])
        parinf.parname = "Iter "+PRING(k+1)+": "+parinf.parname
        whbad = WHERE(parm LT 1.D-4*parinfo[x,y,*].value $
                      AND STRPOS(parinfo[x,y,*].parname,"I(") EQ 0 $ ;; a modifier !!!
                      AND parinfo[x,y,*].fixed EQ 0,Nbad) ;; a modifier !!!
;; Faire un restart selectif...
       IF (Nbad GT 0) THEN parm[whbad] = parinfo[x,y,whbad].value
       fit = MPCURVEFIT( wavobs,Fnuobs0,weights0, $
                          parm,FUNCTION_NAME="specmodel",PARINFO=parinf, $
                          STATUS=status1D,/NODERIVATIVE,FUNCTARGS=extra, $
                          FTOL=ftol,XTOL=1.D-50,GTOL=1.D-50,QUIET=quietFIT, $
                          COVAR=covar,CHI=chi21D,ITMAX=Niter)
      ENDFOR
    ENDIF
    parstrOUT[x,y] = PARSTRUCT(parm,Nbb=Nbb,Nline=Nline,Nband=Nband,NAv=NAv, $
                               Nstar=Nstar)
    
    ;; Error
    sigpar = DBLARR(Nparm)
    FOR i=0,Nparm-1 DO sigpar[i] = SQRT(covar[i,i])

    ;; Final printing
    timin[x,y] = (SYSTIME(1)-t1)/60.
    IF (quiet LE 2) THEN BEGIN
      FOR i=0,Nparm-1 DO IF (parinfo[x,y,i].mpprint EQ 1) THEN $
        PRINT, STRING(parinfo[x,y,i].parname,"(A20)")+" = ", parm[i], " +- ", $
               sigpar[i]
      PRINT, " - status="+PRING(status1D)+", " $
             +PRING((SYSTIME(1)-t0)/60.,F="(F8.1)")+" minutes elapsed."
    ENDIF

  ENDIF ELSE BEGIN
    
    ;; No fit
    sigpar = DBLARR(Nparm)
    status = -1
    chi2 = 0.D

  ENDELSE


  ;; Outputs
  ;;--------
  ;; Final parameters
  errparstrOUT[x,y] = PARSTRUCT(sigpar,Nbb=Nbb,Nline=Nline,Nband=Nband, $
                                NAv=NAv,Nstar=Nstar)
  parinfoOUT = parinfo
  statusOUT[x,y] = status1D
  chi2OUT[x,y] = chi21D

  ;; Save the model
  SPECMODEL, wmod, parm, Fnu_tot1D, FNU_BB=Fnu_BB1D, FNU_STAR=Fnu_star1D, $
             FNU_LINE=Fnu_line1D, FNU_BAND=Fnu_band1D, PABS=Pabs1D, $
             _EXTRA=extra
  IF KEYWORD_SET(modsaveIN) THEN BEGIN
    Fnu_tot[x,y,*] = Fnu_tot1D
    Fnu_BB[x,y,*,*] = Fnu_BB1D
    Fnu_line[x,y,*] = Fnu_line1D
    Fnu_band[x,y,*,*] = Fnu_band1D
    Fnu_star[x,y,*] = Fnu_star1D
    Pabs[x,y,*] = Pabs1D
    fitOUT = {Nw:Nwmod, w:wmod, nu:numod, Fnu_tot:REFORM(Fnu_tot), $
              Fnu_BB:Fnu_BB, Fnu_star:Fnu_star, Fnu_line:Fnu_line, $
              Fnu_band:Fnu_band, Pabs:Pabs}
  ENDIF ELSE fitOUT = 0

  ;; Correct for the absorption of the lines that were fitted first
  FOR i=0,Nline1-1 DO BEGIN
    Pabsline = Pabs1D[CLOSEST(wmod,lines[whline1[i]].wave)]
    parstrOUT[x,y].Iline[whline1[i]] /= Pabsline
    errparstrOUT[x,y].Iline[whline1[i]] /= Pabsline
  ENDFOR

  ;; Plottings
  FOR k=(noplot GE 1),KEYWORD_SET(filepsIN) DO BEGIN
    extra.noplot = 0
    IF (k EQ 1) THEN SETPLOT_PS, /TeX, FILE=fileps[x,y]
      SPECMODEL, wmod, parm, Fnu_tot1D, FNU_BB=Fnu_BB1D, FNU_STAR=Fnu_star1D, $
                 FNU_LINE=Fnu_line1D, FNU_BAND=Fnu_band1D, PABS=Pabs1D, $
                 _EXTRA=extra
    IF (k EQ 1) THEN $
      ENDPS, FILE=fileps[x,y], PDF=KEYWORD_SET(pdfIN), EPS=KEYWORD_SET(epsIN), $
             NOX=KEYWORD_SET(noX) 
  ENDFOR

ENDIF
ENDFOR

  ;; Savings (once every column)
  IF (KEYWORD_SET(filexdrIN)) THEN BEGIN
    parstr = REFORM(parstrOUT)
    errparstr = REFORM(errparstrOUT)
    fit = fitOUT
    status = REFORM(statusOUT)
    chi2 = REFORM(chi2OUT)
    SAVE, FILE=filexdrIN, Nx, Ny, Nobs, wavobs, nuobs, Fnuobs, dFnuobs, $
                          weights, mask, parstr, errparstr, parinfo, $
                          status, chi2, fit, instrument, timin, cont, lines, $
                          bands, extinct
    PRINT, " - File "+filexdrIN+" has been written."
  ENDIF 

ENDFOR
;; BIG LOOP

IF (quiet LE 1) THEN PRINT, STRJOIN(REPLICATE("-",80))


END
