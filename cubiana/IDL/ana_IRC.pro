src = 'M82'
data_path = '/Users/dhu/Data/PAHPedia/'+src+'/'
fit_path = data_path+'fit/'

Fnu = MRDFITS(data_path+src+"_IRC_0.fits",0,hdr)
RESTORE, fit_path+"fit_IRC_0.xdr"

Nmc = 20

for i =0,Nmc do begin
	;; Sum up sub-components
	RESTORE, fit_path+"fit_IRC_"+PRING(i)+".xdr" ;; loads the fit results
	parstr.Iline *= !MKS.Jy ;; conversion [MJy.Hz/sr] --> [MW/m2/sr]
	parstr.Iband *= !MKS.Jy ;; conversion [MJy.Hz/sr] --> [MW/m2/sr]
	I33 = parstr.Iband[WHERE(bands.label EQ "Main 3.3")]
	I34 = parstr.Iband[WHERE(bands.label EQ "Main 3.4")]
	
	;; Header
	hdr0 = hdr
	SXDELPAR, hdr0, "NAXIS3"
	SXDELPAR, hdr0, "COMMENT"
	SXADDPAR, hdr0, "COMMENT", "Fitted band intensity, in MW/m2/sr"

	;; Write fits files of the main features
	MWRFITS, REFORM(I33), fit_path+"Intensities/I33/I33_"+PRING(i)+".fits", hdr0, /CREATE
	MWRFITS, REFORM(I34), fit_path+"Intensities/I34/I34_"+PRING(i)+".fits", hdr0, /CREATE
endfor

end
