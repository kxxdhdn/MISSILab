from astropy.io import fits
from astropy.nddata import CCDData, Cutout2D
from astropy import units as u

filename = 'dh_SL2_cube'
with fits.open(filename+'.fits') as hdul0:
	#data0 = hdul0[0].data
	hdr0 = hdul0[0].header
	NAXIS1 = hdr0['NAXIS1']
	NAXIS2 = hdr0['NAXIS2']
	NAXIS3 = hdr0['NAXIS3']
	CRPIX1 = hdr0['CRPIX1']
	CRPIX2 = hdr0['CRPIX2']
xl = 2 # set the number to cut (x left side)
xr = 2
yl = 0
yr = 0
nx = NAXIS1 - xl - xr
ny = NAXIS2 - yl - yr
xc = int(nx / 2) # NO "+1" because of 0-indexed 
yc = int(ny / 2)
position = (xc, yc) # (x, y)
size = (ny, nx) # (ny, nx) centering convention: small x & y have the priority when even number
ccd = CCDData.read(filename+'.fits')
print("---------- Ignore the info above ----------")
print("nx, ny, xc, yc = ", nx, ny, xc, yc)
#cutout = ccd[:,:ny,:nx]
cutout = ccd[:,yl:ny+yl,xl:nx+xl]
for i in range(NAXIS3):
	cutout[i,:,:].data = Cutout2D(ccd[i,:,:], position, size)
print(ccd.shape)
print(cutout.shape)
hdr0['CRPIX1'] -= xl
hdr0['CRPIX2'] -= yl
cutout.meta = hdr0
cutout.write(filename+'_c.fits', overwrite=True)
