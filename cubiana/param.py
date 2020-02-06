import os

src = 'M82'
Nmc = 60

## Current dir
##-------------
path_cur = os.getcwd()+'/'
path_idl = path_cur+'/IDL/'

## Root dir
##----------
root = '/Users/dhu/Data/'

## Data dir
##----------
path_irs = root+'Spitzer/data/'+src+'/'
path_irc = root+'AKARI/data/'+src+'/'
path_phot = root+'Photometry/'+src+'/'
path_ker = root+'Kernels/'

## IRS data (CUBISM)
##-------------------
fits_sl2 = []
fits_sl3 = []
fits_sl1 = []
fits_ll2 = []
fits_ll3 = []
fits_ll1 = []
fits_irs = []
chnl = ['SL2', 'SL3', 'SL1', 'LL2', 'LL3', 'LL1']
lab_sl = ['_04', '_06', '_08', '_09']
lab_ll = ['_04', '_05', '_06', '_08', '_09']
# lab_sl = ['']
# lab_ll = ['']
for i, ch in enumerate(chnl):
	## SL
	if i//3==0:
		for t in lab_sl:
			f = path_irs+src+t+'_'+ch
			if ch=='SL2':
				fits_sl2.append(f)
			if ch=='SL3':
				fits_sl3.append(f)
			if ch=='SL1':
				fits_sl1.append(f)
	## LL
	else:
		for t in lab_ll:
			f = path_irs+src+t+'_'+ch
			if ch=='LL2':
				fits_ll2.append(f)
			if ch=='LL3':
				fits_ll3.append(f)
			if ch=='LL1':
				fits_ll1.append(f)
fits_irs.append(fits_sl2)
fits_irs.append(fits_sl3)
fits_irs.append(fits_sl1)
fits_irs.append(fits_ll2)
fits_irs.append(fits_ll3)
fits_irs.append(fits_ll1)

## PSF
##-----
fits_ker = []
psf = [2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6.]
# psf_ref = 'IRAC_5.8' # 2.11 (< LL1)
# psf_ref = 'IRAC_8.0'# 2.82 (< LL1)
psf_ref = 'Gauss_06.0'
# psf_ref = 'MIPS_24' # 6.43
# psf_ref = 'WISE_MAP_11.6' # 6.60
# psf_ref = 'WISE_MAP_22.1' # 11.89
for p in psf:
	fits_ker.append(path_ker+'Kernel_HiRes_Gauss_0'+\
		str(p)+'_to_'+psf_ref)

## Outputs
##---------
path_out = root+'PAHpedia/'+src+'/'
## See also IDL/conv_prog.pro
path_tmp = root+'PAHpedia/tmp/'
if not os.path.exists(path_tmp):
	os.makedirs(path_tmp)
## See also IDL/convolve_image.pro
path_slices = path_tmp+'slices/'
if not os.path.exists(path_slices):
	os.makedirs(path_slices)
## See also IDL/conv_prog.pro
csv_ker = path_tmp+'kernels'
## IRC cube building
path_build = path_tmp+'cubuild/'
if not os.path.exists(path_build):
	os.makedirs(path_build)

## IRC data
##----------
slits = ['Ns', 'Nh']
obsid = [
	['3390001.1','F011100297_N002'], 
	['3390001.2','F011100338_N002'], 
	['3390002.1','F011100795_N002'], 
	['3390002.2','F011176073_N002'], 
	['3390003.1','F011100379_N002'], 
]
fits_irc = []
parobs = []
for obs in obsid:
	for s in slits:
		fits_irc.append(path_build+obs[0]+'_'+s)
		parobs.append([obs[0], s, obs[1]])

## Calibrations
##--------------
phot = 'IRAC3' # syn_phot ref
phot0 = 'IRAC3_SINGS' # phot ref
path_cal = path_out+'calib/'
if not os.path.exists(path_cal):
	os.makedirs(path_cal)
fits_phot = path_phot+src+'_cal_'+phot
fits_phot0 = path_phot+src+'_cal_'+phot0

## Tests
##-------
path_tests = root+'PAHpedia/tests/'

verbose = False
