from astropy.io import fits
from astropy import units as u

xl = 2 # set cut pixel number (x left border)
xr = 2
yl = 0
yr = 0
filename = '0dh_SL2_cube'
tablename = 'wvl'
with fits.open(filename+'.fits') as hdul0:
	data0 = hdul0[0].data
	hdr0 = hdul0[0].header
	NAXIS1 = hdr0['NAXIS1']
	NAXIS2 = hdr0['NAXIS2']
	NAXIS3 = hdr0['NAXIS3']
	hdr0['CRPIX1'] -= xl # ref ax
	hdr0['CRPIX2'] -= yl # ref ay
	# ref ra & dec not changed
	table = hdul0[1].data # table = [array([[*],[*], ...],)] list/array/list/list
	wvl = []
	for i in range(NAXIS3):
		wvl.append(table[0][0][i][0])
	hdr1 = hdul0[1].header
primary_hdu = fits.PrimaryHDU(header=hdr0, data=data0[:,yl:NAXIS2-yr,xl:NAXIS1-xr])
hdul = fits.HDUList(primary_hdu)
hdu = fits.ImageHDU(data=wvl, header=hdr1, name=tablename) # add table
hdul.append(hdu)
hdul.writeto(filename+'_c.fits', output_verify='ignore', overwrite=True)
