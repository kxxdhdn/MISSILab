from astropy.io import fits
from astropy import units as u
from astropy import wcs
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors


def mapview(x, y, wvl, data, radius):
	"""
	show the map and choose 2 points to show spectra

	--- INPUT ---
	wvl      
	data     
	radius   int type
	"""
	wvlnum = 10
	fig = plt.figure("Choose 2 points to compare their spectra")
	cnorm = None # Linear
#	cnorm = colors.PowerNorm(gamma=0.5) # Square root
#	cnorm = colors.SymLogNorm(linthresh=1., linscale=1., clip=False) # Logarithmic
	ax = plt.subplot()
	pcm = ax.pcolormesh(x, y, data[wvlnum,:,:], norm=cnorm, cmap='rainbow')
	fig.colorbar(pcm, extend='max')
	ax.set_title("IRS map at {} micron".format(wvl[wvlnum]))
	ax.set_xlabel("x")
	ax.set_ylabel("y")
	# left click to choose 2 points
	pts = plt.ginput(n=2, show_clicks=True)
	(x1, y1) = np.array(pts[0]).astype(int)
	(x2, y2) = np.array(pts[1]).astype(int)
	# average flux in a circle of radius size
	print("radius = ", radius)
	data1 = data[:,y1,x1]
	data2 = data[:,y2,x2]
	n = 1
	for dx in range(radius+1):
		for dy in range(radius+1):
			p2 = dx*dx+dy*dy
			r2 = radius*radius
			if (p2!=0 and p2<=r2):
				#print("(dx, dy) = ", dx, dy)
				data1 += data[:,y1-dy,x1-dx] + data[:,y1-dy,x1+dx] + data[:,y1+dy,x1-dx] + data[:,y1+dy,x1+dx]
				data2 += data[:,y2-dy,x2-dx] + data[:,y2-dy,x2+dx] + data[:,y2+dy,x2-dx] + data[:,y2+dy,x2+dx]
				n += 4
	print("flux mean of pixels centered at ({0}, {1}) & ({2}, {3})".format(x1, y1, x2, y2))
	# plot spectrum at chosen points
	fig2 = plt.figure("spectra")
	ax2 = plt.subplot()
	ax2.plot(wvl, data1/n, 'c-', linewidth=.5, label='1st point')
	ax2.plot(wvl, data2/n, 'm-', linewidth=.5, label='2nd point')
	ax2.set_xlabel('Wavelength (micron)')
	ax2.set_ylabel(r'$F_{\nu} (MJy)$')
#	ax2.set_yscale('log', nonposy='clip')
	ax2.legend(loc='best')
	
	fig2.savefig('spectra.png')
	plt.show()

def bimapview(x, y, wvl, data, data0, radius):
	"""
	show the map and choose 1 points to show spectra

	--- INPUT ---
	wvl      
	data     
	radius   int type
	"""
	wvlnum = 10
	fig1 = plt.figure("Choose 2 points to compare their spectra", figsize=(14,5))
	fig1.suptitle("IRS maps at {} micron".format(wvl[wvlnum]))
	cnorm = None # Linear
#	cnorm = colors.PowerNorm(gamma=0.5) # Square root
#	cnorm = colors.SymLogNorm(linthresh=1., linscale=1., clip=False) # Logarithmic
	ax11 = plt.subplot(121)
	pcm = ax11.pcolormesh(x, y, data[wvlnum,:,:], norm=cnorm, cmap='rainbow')
	fig1.colorbar(pcm, ax=ax11, extend='max')
	ax11.set_title("dh")
	ax11.set_xlabel("x")
	ax11.set_ylabel("y")
	ax12 = plt.subplot(122)
	pcm = ax12.pcolormesh(x, y, data0[wvlnum,:,:], norm=cnorm, cmap='rainbow')
	fig1.colorbar(pcm, ax=ax12, extend='max')
	ax12.set_title("sh")
	ax12.set_xlabel("x")
	ax12.set_ylabel("y")
	# left click to choose 1 points
	pts = plt.ginput(n=1, show_clicks=True)
	(x1, y1) = np.array(pts[0]).astype(int)
	# average flux in a circle of radius size
	print("radius = ", radius)
	data1 = data[:,y1,x1]
	data2 = data0[:,y1,x1]
	n = 1
	for dx in range(radius+1):
		for dy in range(radius+1):
			p2 = dx*dx+dy*dy
			r2 = radius*radius
			if (p2!=0 and p2<=r2):
#				print("(dx, dy) = ", dx, dy)
				data1 += data[:,y1-dy,x1-dx] + data[:,y1-dy,x1+dx] + data[:,y1+dy,x1-dx] + data[:,y1+dy,x1+dx]
				data2 += data0[:,y1-dy,x1-dx] + data[:,y1-dy,x1+dx] + data[:,y1+dy,x1-dx] + data[:,y1+dy,x1+dx]
				n += 4
	print("flux mean of pixels centered at ({0}, {1})".format(x1, y1))
	# plot spectrum at chosen points
	fig2 = plt.figure("spectra")
	ax2 = plt.subplot()
	ax2.plot(wvl, data1/n, 'c-', linewidth=.5, label='dh')
	ax2.plot(wvl, data2/n, 'm-', linewidth=.5, label='sh')
	ax2.set_xlabel('Wavelength (micron)')
	ax2.set_ylabel(r'$F_{\nu} (MJy)$')
#	ax2.set_yscale('log', nonposy='clip')
	ax2.legend(loc='best')
	
	fig2.savefig('spectra2.png')
	plt.show()


filename = 'SL2_cube_o'
#filename0 = 'SL2_cube_b'
filename0 = 'sh_SL2_cube'
with fits.open(filename+'.fits') as hdul:
	data = hdul[0].data
	hdr = hdul[0].header
	NAXIS1 = hdr['NAXIS1']
	NAXIS2 = hdr['NAXIS2']
	NAXIS3 = hdr['NAXIS3']
	wvl = hdul[1].data[0][0][:,0]
#	wvl = hdul[1].data# for rewitten header
with fits.open(filename0+'.fits') as hdul0:
	data0 = hdul0[0].data

x = np.arange(hdr['NAXIS1'])
y = np.arange(hdr['NAXIS2'])
i = 0
for i in range(10):
	i += 1
#	mapview(x, y, wvl, data, 2)
	bimapview(x, y, wvl, data, data0, 1)
