import os

verbose = False

src = 'M82'
Nmc = 0

## Current dir
##-------------
path_cur = os.getcwd()+'/'
path_par = os.path.dirname(os.path.abspath(__file__))+'/' # param file path

## IDL dir
##---------
path_idl = path_par+'idlib/'

## Root dir
##----------
path_root = '/Users/dhu/Data/'
## Data dir
path_irc = path_root+'AKARI/data/'+src+'/'
path_phot = path_root+'Photometry/'+src+'/'
path_ker = path_root+'Kernels/'

## Convolution
##-------------
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

## Tmp files
##-----------
path_tmp = path_root+'PAHPedia/tmp/'
if not os.path.exists(path_tmp):
    os.makedirs(path_tmp)
csv_ker = path_tmp+'kernelist' # idlib/conv_prog.pro

## Outputs
##---------
path_out = path_root+'PAHPedia/tests/Y12spec/'

## IRC cube building
path_build = path_tmp+'cubuild/'
if not os.path.exists(path_build):
    os.makedirs(path_build)

## IRC data
##----------
slits = ['Ns', 'Nh']
obsid = [
    ['3390001.1','F011100297_N002'], 
    # ['3390001.2','F011100338_N002'], 
    ['3390002.1','F011100795_N002'], 
    # ['3390002.2','F011176073_N002'], 
    ['3390003.1','F011100379_N002'], 
    # ['5125401.1', ]
    # ['5125402.1', ]
    # ['5125403.1', ]
    # ['5125404.1', ]
    # ['5125405.1', ]
    # ['5125406.1', ]
]
fits_irc = []
parobs = []
for obs in obsid:
    for s in slits:
        fits_irc.append(path_build+obs[0]+'_'+s)
        parobs.append([obs[0], s, obs[1]])

## Calibrations
##--------------
phot = 'IRAC1' # syn_phot ref
path_cal = path_out+'calib/'
if not os.path.exists(path_cal):
    os.makedirs(path_cal)

## Tests
##-------
path_tests = path_root+'PAHPedia/tests/'
