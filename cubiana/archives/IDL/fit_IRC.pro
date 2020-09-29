;; Path of the start-up file
; @$idlib/idl_init

src = 'M82'
z = 0.00068
data_path = '/Users/dhu/Data/PAHPedia/'+src+'/'
fit_path = data_path+'fit/'

Nmc  = 20

for i=0,Nmc do begin
	;; Read one spectral cube. Data must be in Jy/pixel
	;; Read one spectral cube. Data is in MJy/sr
	Fnu = MRDFITS(data_path+src+"_IRC_"+pring(i)+".fits",0,hdr)
	wave = MRDFITS(data_path+src+"_IRC_"+pring(i)+".fits",1)
	dFnu = MRDFITS(data_path+src+"_IRC_unc.fits",0,hdr)
	
	;; Restrain the wavelength range
	wh = WHERE(wave GE 2.8 AND wave LE 21.,Nw)
	Fnu = Fnu[*,*,wh]
	dFnu = dFnu[*,*,wh]
	wave = wave[wh]/(1+z)
	siz = SIZE(Fnu)
    Nx = siz[1]
    Ny = siz[2]
    mask = BYTARR(Nx,Ny)
    FOR j=0,Nw-1 DO BEGIN
      wh = WHERE(FINITE(Fnu[*,*,j]) EQ 0 OR dFnu[*,*,j] LE 0,N)
      IF (N GT 0) THEN mask[wh] = 1 
    ENDFOR
    wh = WHERE(dFnu/Fnu LT 1.D-3,N)
    IF (N GT 0) THEN dFnu[wh] = Fnu[wh]*1.D-3
    wh = WHERE(FINITE(dFnu) EQ 0,N)
    IF (N GT 0) THEN dFnu[wh] = Fnu[wh]*1.D-1

	;; Fit settings
	lines = ['Bra'] ;; see Table 1, page 2 of notes.pdf
	Nband = 2+31 ;; see Table 3, page 11 of notes.pdf
	ijig = [0,1]
	fixband = REPLICATE(1,Nband)
	fixband[ijig] = 0
	Ncont = 2 ;; number of black body components
	Tini = [180.D,50.D] ;; initial temperature [K]
	bbtype = ["crb","crb"] ;; amorphous carbon and silicates
	
	;; Actual fit
        IF (i EQ 0) THEN fileps = fit_path+"Figures/fit_IRC.eps" ELSE fileps = 0
	MILES, wave, Fnu, dFnu, MASK=mask, INSTRUMENT="Ns-SL", BBTYPE=bbtype, TEMPBB=Tini, /BANDS, $
	       LINES=lines, /OLD_STARS, FIXCENTERBAND=fixband, FIXSIGMABAND=fixband, $
	       NOPLOT=2, QUIET=2, SAVE=fit_path+"fit_IRC_"+PRING(i)+".xdr", $
	       COPY=fileps, /PDF

endfor

end
