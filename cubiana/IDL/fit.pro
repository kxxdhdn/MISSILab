;; Path of the start-up file
;@/Users/dhu/PhD/pahpedia/SWING/idl_init

source = "M82_O"
z = 0.00068
data_path = "../data/"+source+"/output/"
fit_path = "../fitting/"+source+"/"

Nmc  = 5

for i=0,Nmc do begin
	;; Read one spectral cube. Data must be in Jy/pixel
	;; Read one spectral cube. Data is in MJy/sr
	Fnu = MRDFITS(data_path+source+"_"+pring(i)+".fits",0,hdr)
	wave = MRDFITS(data_path+source+"_"+pring(i)+".fits",1)
	dFnu = MRDFITS(data_path+source+"_unc.fits",0,hdr)
	
	;; Restrain the wavelength range
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
	IF (i EQ 0) THEN fileps = fit_path+"Figures/fit_cube.eps" ELSE fileps = 0
	MILES, wave, Fnu, dFnu, INSTRUMENT="SL", BBTYPE=bbtype, TEMPBB=Tini, /BANDS, $
	       LINES=lines, /OLD_STARS, FIXCENTERBAND=fixband, FIXSIGMABAND=fixband, $
	       NOPLOT=2, QUIET=2, SAVE=fit_path+"fit_cube_"+PRING(i)+".xdr", $
	       COPY=fileps, /PDF

endfor

end
