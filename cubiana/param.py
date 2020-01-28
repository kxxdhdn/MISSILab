import os

src = 'M82'
Nmc = 2

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
fits_irs = []
fits_unc_irs = []
chnl = ['SL2', 'SL1', 'LL2']#, 'LL1']
Tobs = ['04', '06', '08', '09']
# Tobs = None
for ch in chnl:
	if Tobs is not None:
		for t in Tobs:
			fits_irs.append(path_irs+src+'_'+t+'_'+ch)
			fits_unc_irs.append(path_irs+src+'_'+t+'_'+ch+'_unc')
	else:
		fits_irs.append(path_irs+src+'_'+ch)
		fits_unc_irs.append(path_irs+src+'_'+ch+'_unc')

## IRC data
##----------


## PSF
##-----
fits_ker = []
psf = [2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6.]
# psf_ref = 'IRAC_5.8' # 2.11 (< LL1)
# psf_ref = 'IRAC_8.0'# 2.82 (< LL1)
# psf_ref = 'Gauss_06.0'
psf_ref = 'MIPS_24' # 6.43
# psf_ref = 'WISE_MAP_11.6' # 6.60
# psf_ref = 'WISE_MAP_22.1' # 11.89
for p in psf:
	fits_ker.append(path_ker+'Kernel_HiRes_Gauss_0'+\
		str(p)+'_to_'+psf_ref)

## Outputs
##---------
path_out = root+'PAHpedia/'+src+'/'

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
path_test = root+'PAHpedia/tests/'
path_tmp = root+'PAHpedia/tmp/'
if not os.path.exists(path_tmp):
	os.makedirs(path_tmp)
verbose = False
