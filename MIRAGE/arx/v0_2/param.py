import os

verbose = False

src = 'M82'
Nmc = 100

## Current dir
##-------------
path_cur = os.getcwd()+'/'
path_par = os.path.dirname(os.path.abspath(__file__))+'/' # param file path

## IDL dir
##---------
path_idl = path_par+'idlib/'

## Root dir
##----------
root = '/Users/dhu/Data/'

## Data dir
##----------
path_irs = root+'Spitzer/'+src+'/'
path_irc = root+'AKARI/'+src+'/'
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
chnl = ['SL2','SL1','LL2','LL1']
# chnl = ['SL2', 'SL3', 'SL1', 'LL2', 'LL3', 'LL1']
Nch = len(chnl) # Number of chnl (ordered as above) used
lab_sl = ['04','06S','06N','08','08c','09N3','09N2']
lab_ll = ['03','04','05','06','08','09N3','09N5','09N2']
# lab_sl = ['']
# lab_ll = ['']
for ch in chnl:
    ## SL
    if ch[:2]=='SL':
        for t in lab_sl:
            f = path_irs+src+'_'+t+'_'+ch
            if ch=='SL2':
                fits_sl2.append(f)
            if ch=='SL3':
                fits_sl3.append(f)
            if ch=='SL1':
                fits_sl1.append(f)
    ## LL
    elif ch[:2]=='LL':
        for t in lab_ll:
            f = path_irs+src+'_'+t+'_'+ch
            if ch=='LL2':
                fits_ll2.append(f)
            if ch=='LL3':
                fits_ll3.append(f)
            if ch=='LL1':
                fits_ll1.append(f)

if 'SL2' in chnl:
    fits_irs.append(fits_sl2)
if 'SL3' in chnl:
    fits_irs.append(fits_sl3)
if 'SL1' in chnl:
    fits_irs.append(fits_sl1)
if 'LL2' in chnl:
    fits_irs.append(fits_ll2)
if 'LL3' in chnl:
    fits_irs.append(fits_ll3)
if 'LL1' in chnl:
    fits_irs.append(fits_ll1)

## PSF
##-----
fits_ker = []
psf = [1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6.]
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
path_out = root+'PAHPedia/'+src+'/'
## See also IDL/conv_prog.pro
path_tmp = root+'PAHPedia/tmp/'
if not os.path.exists(path_tmp):
	os.makedirs(path_tmp)
## See also IDL/convolve_image.pro
path_slices = path_tmp+'slices/'
if not os.path.exists(path_slices):
	os.makedirs(path_slices)
## See also IDL/conv_prog.pro
csv_ker = path_tmp+'kernelist'
## IRC cube building
path_build = path_tmp+'cubuild/'
if not os.path.exists(path_build):
	os.makedirs(path_build)

## IRC data
##----------
slits = ['Ns', 'Nh']
obsid = [
    ['3390001.1','F011100297_N002','NG'], # A1-6, H1-2
    ['3390001.2','F011100338_N002','NG'], # 
    ['3390002.1','F011100795_N002','NG'], # C1-6, I1-2
    # ['3390002.2','F011176073_N002','NG'], # Matching failed
    ['3390003.1','F011100379_N002','NG'], # B1-6, G1-2
    ['5124077.1','F007174142_N002','NG'], # cap
    ['5125401.1','F010117172_N002','NG'], # J1-2, E1-3
    # ['5125402.1','F010117172_N002','NP'], # Matching failed
    ['5125403.1','F010116924_N002','NG'], # 
    # ['5125404.1','F010117338_N002','NP'], # Matching failed
    ['5125405.1','F010116950_N002','NG'], # D1-3, F1-2
    # ['5125406.1','F010117086_N002','NP'], # NP has 68 wvl instead of 259
]
fits_irc = []
parobs = []
for obs in obsid:
	for s in slits:
		fits_irc.append(path_build+obs[0]+'_'+s)
		parobs.append([obs[0], s, obs[1]])

## Calibrations
##--------------
phot = 'IRAC4' # syn_phot ref
phot0 = 'IRAC4_SINGS' # phot ref
path_cal = path_out+'calib/'
if not os.path.exists(path_cal):
	os.makedirs(path_cal)
fits_phot = path_phot+src+'_'+phot
fits_phot0 = path_phot+src+'_'+phot0

## Tests
##-------
path_tests = root+'PAHPedia/tests/'

verbose = False
