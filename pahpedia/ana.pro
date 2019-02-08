Fnu = MRDFITS("../../data/m82/reprojection/m82_reproj_0.fits",0,hdr)
RESTORE, "fit_cube_0.xdr"

Nmc = 10
cube62 = DBLARR(Nx,Ny,Nmc)
cube77 = DBLARR(Nx,Ny,Nmc)
cube86 = DBLARR(Nx,Ny,Nmc)
cube112 = DBLARR(Nx,Ny,Nmc)
cube127 = DBLARR(Nx,Ny,Nmc)
cube170 = DBLARR(Nx,Ny,Nmc)

for i =0,Nmc do begin
	;; Sum up sub-components
	RESTORE, "fit_cube_"+PRING(i)+".xdr" ;; loads the fit results
	parstr.Iline *= !MKS.Jy ;; conversion [Jy.Hz] --> [W/m2]
	parstr.Iband *= !MKS.Jy ;; conversion [Jy.Hz] --> [W/m2]
	I62 = parstr.Iband[WHERE(bands.label EQ "Main 6.2 (1)")] $
		+ parstr.Iband[WHERE(bands.label EQ "Main 6.2 (2)")]
	I77 = parstr.Iband[WHERE(bands.label EQ "Main 7.7 (1)")] $
		+ parstr.Iband[WHERE(bands.label EQ "Main 7.7 (2)")] $
		+ parstr.Iband[WHERE(bands.label EQ "Plateau 7.7")]
	I86 = parstr.Iband[WHERE(bands.label EQ "Main 8.6")]
	I112 = parstr.Iband[WHERE(bands.label EQ "Main 11.2")] $
		+ parstr.Iband[WHERE(bands.label EQ "Plateau 11.3")]
	I127 = parstr.Iband[WHERE(bands.label EQ "Main 12.7 (1)")] $
		+ parstr.Iband[WHERE(bands.label EQ "Main 12.7 (2)")]
	I170 = parstr.Iband[WHERE(bands.label EQ "Plateau 17.0")] $
		+ parstr.Iband[WHERE(bands.label EQ "Small 16.4")] $
		+ parstr.Iband[WHERE(bands.label EQ "Small 17.4")] $
		+ parstr.Iband[WHERE(bands.label EQ "Small 17.8")]
	
	;; Header
	hdr0 = hdr
	SXDELPAR, hdr0, "NAXIS3"
	SXDELPAR, hdr0, "COMMENT"
	SXADDPAR, hdr0, "COMMENT", "Fitted band intensity, in MW/m2/sr"

	;; Write fits files of the main features
	MWRFITS, REFORM(I62), "Intensity/I62/I62_"+PRING(i)+".fits", hdr0, /CREATE
	MWRFITS, REFORM(I77), "Intensity/I77/I77_"+PRING(i)+".fits", hdr0, /CREATE
	MWRFITS, REFORM(I86), "Intensity/I86/I86_"+PRING(i)+".fits", hdr0, /CREATE
	MWRFITS, REFORM(I112), "Intensity/I112/I112_"+PRING(i)+".fits", hdr0, /CREATE
	MWRFITS, REFORM(I127), "Intensity/I127/I127_"+PRING(i)+".fits", hdr0, /CREATE
	MWRFITS, REFORM(I170), "Intensity/I170/I170_"+PRING(i)+".fits", hdr0, /CREATE
endfor

end
