SUGGESTS

In the "../scripts" path enter the command "python reproj.py" and then follow the instructions.

Coordinates of center galaxy: (204 15 10, -29 -51 -52)

Coordinates for zooming into the center:
xmin = 204 13 50
xmax = 204 16 17
ymin = -29 -53 -10
ymax = -29 -50 -37

It is possible that you should install extra python packages like "reproject" (e.g. with Anaconda you can simply use "conda install reproject" to make it.)



SCRIPTS

	PATH: 
	../scripts

	CONTENT: 
	hdupy.py	Header Data Unit operations
	visupy.py	spatial and spectral visualization of fits file
	test.py		general package ("py" as postfix) test
	unc.py		uncertainty estimation of the background sky by comparing with the *_unc.fits data.
	dec.py		separate the raw fits files to fits files with data of a single wavelength.
	conv.py		call IDL to convolve images with PSF kernels
	reproj.py	reproject galaxy observations so that we can recombine images of all wavelengths to a single fits file and later plot the integrated spectrum. (e.g. using visualisation.py)



DATA

M83 raw data: 
4 fits files (3D) of tunnels sl2, sl1, ll1, ll2 and their corresponding uncertainty fits files (2D) 
PATH: 
../data/m83/
FILENAMES: 
m83_sl2_cube.fits

M83 single wavelength images:
369 fits files (2D)
PATH:
../data/m83/slices/
FILENAMES: 
m83_sl2_0.fits

M83 convolved images:
369 fits files (2D)
PATH: 
../data/m83/convolved/
FILENAMES: 
m83_sl2_0_conv.fits

PSF kernels: 
20 fits files (2D)
PATH: 
../data/kernels/
FILENAMES: 
Kernel_HiRes_Gauss_01.5_to Gauss_10.0fits