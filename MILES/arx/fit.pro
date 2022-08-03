;; Path of the start-up file
;@/Users/dhu/stage/fitting/SWING/idl_init
Nmc  = 20

for i=0,Nmc do begin
	;; Read one spectral cube. Data must be in Jy/pixel
;	Fnu = MRDFITS(strcompress("../../data/m82/reprojection/m82_reproj_"+string(boucle)+".fits",/remove),0,hdr)
;	wave = MRDFITS(strcompress("../../data/m82/reprojection/m82_reproj_"+string(boucle)+".fits",/remove),1)
;	dFnu = MRDFITS("m82_reproj_unc.fits",0,hdr)
	Fnu = MRDFITS("../../data/m82/reprojection/m82_reproj_"+pring(i)+".fits",0,hdr)
	wave = MRDFITS("../../data/m82/reprojection/m82_reproj_"+pring(i)+".fits",1)
	dFnu = MRDFITS("../../data/m82/reprojection/m82_reproj_unc.fits",0,hdr)
	
	;; Restrain the wavelength range
	z = 0.000677
	wh = WHERE(wave LE 21.)
	Fnu = Fnu[*,*,wh]
	dFnu = dFnu[*,*,wh]
	wave = wave[wh]/(1+z)
	
	;; Fit settings
	lines = ["H2S7","ArII","H2S5","ArIII1","H2S3","SIV","H2S2", $
	         "NeII","NeIII1","H2S1","SIII1"] ;; see Table 1, page 2 of notes.pdf
	Nband = 2+31 ;; see Table 3, page 11 of notes.pdf
	ijig = 2 + [4,10,13,17]
	fixband = REPLICATE(1,Nband)
	fixband[ijig] = 0
	Ncont = 4 ;; number of black body components
	Tini = [80.D,150.D,180.D,50.D] ;; initial temperature [K]
	bbtype = ["crb","sil_norm","crb","crb"] ;; amorphous carbon and silicates
	
	;; Actual fit
	IF (i EQ 0) THEN fileps = "Figures/fit_cube.eps" ELSE fileps = 0
	MILES, wave, Fnu, dFnu, INSTRUMENT="SL", BBTYPE=bbtype, TEMPBB=Tini, /BANDS, $
	       LINES=lines, /OLD_STARS, FIXCENTERBAND=fixband, FIXSIGMABAND=fixband, $
	       NOPLOT=2, QUIET=2, SAVE="fit_cube_"+PRING(i)+".xdr", $
	       COPY=fileps, /PDF
	
	;;; Sum up sub-components
	;RESTORE, "fit_cube_0.xdr" ;; loads the fit results
	;parstr.Iline *= !MKS.Jy ;; conversion [Jy.Hz] --> [W/m2]
	;parstr.Iband *= !MKS.Jy ;; conversion [Jy.Hz] --> [W/m2]
	;I62 = parstr.Iband[WHERE(bands.label EQ "Main 6.2 (1)")] $
	;    + parstr.Iband[WHERE(bands.label EQ "Main 6.2 (2)")]
	;I77 = parstr.Iband[WHERE(bands.label EQ "Main 7.7 (1)")] $
	;    + parstr.Iband[WHERE(bands.label EQ "Main 7.7 (2)")] $
	;    + parstr.Iband[WHERE(bands.label EQ "Plateau7 .7")]
	;I86 = parstr.Iband[WHERE(bands.label EQ "Main 8.6")]
	;I112 = parstr.Iband[WHERE(bands.label EQ "Main 11.2")] $
	;     + parstr.Iband[WHERE(bands.label EQ "Plateau 11.3")]
	;I127 = parstr.Iband[WHERE(bands.label EQ "Main 12.7 (1)")] $
	;     + parstr.Iband[WHERE(bands.label EQ "Main 12.7 (2)")]
	;
	;;; Header
	;hdr0 = hdr
	;SXDELPAR, hdr0, "NAXIS3"
	;SXDELPAR, hdr0, "COMMENT"
	;SXADDPAR, hdr0, "COMMENT", "Fitted band intensity, in MW/m2/sr"
	;
	;;; Write fits files of the main features
	;MWRFITS, REFORM(I62), "I62_"+pring(i)+".fits", hdr0, /CREATE
	;MWRFITS, REFORM(I77), "I77_"+pring(i)+".fits", hdr0, /CREATE
	;MWRFITS, REFORM(I86), "I86_"+pring(i)+".fits", hdr0, /CREATE
	;MWRFITS, REFORM(I112), "I112_"+pring(i)+".fits", hdr0, /CREATE
	;MWRFITS, REFORM(I127), "I127_"+pring(i)+".fits", hdr0, /CREATE
endfor

end
