#import sys
#sys.path.append('..')
#print(sys.path)
import numpy as np
import matplotlib.pyplot as plt
import csv
from os import system
from myfunc import kern


plot_save_path = '../plots/m83/'
chnl = ['sl2', 'sl1', 'll2', 'll1']
#chnl = ['akari']
Nc = np.size(chnl)
psf = [01.0, 01.5, 02.0, 02.5, 03.0, 03.5, 04.0, 04.5, 05.0]
#, 05.5, 06.0, 06.5, 07.0, 07.5, 08.0, 08.5, 09.0, 09.5, 10.0, 11.0, 12.0]


def choker(index, wvli, wvl):
	k = 0
	choices = []
	for i in range(Nc):
		sigma_lam = kern(wvli[i])[0]
		##choose kernels
		for j in range(np.size(wvli[i])):
			if i == Nc: #attention! not adapted to all (Nc-1 with akari, in decomp.py)
				kern_name = 'Kernel_HiRes_Gauss_0' + str(psf[8]) + '_to_Gauss_06.0'
			else:
				for p in range(np.size(psf)):
					dp = sigma_lam[j] - psf[p]
					dp1 = sigma_lam[j] - psf[p+1]
					if (dp > 0 and dp1 > 0):
						flag = 0
					elif (dp > 0 and dp1 < 0):
						if dp + dp1 <= 0:
							flag = 1
							break
						else:
							p += 1
							flag = 1
							break
					else:
						print('Maybe need smaller PSF?')
						flag = 1
						break
				if (flag==0):
					print('Maybe need bigger PSF? (In this case IndexError will appear first at psf[p+1])')
				kern_name = 'Kernel_HiRes_Gauss_0' + str(psf[p]) + '_to_Gauss_06.0'
			im_name = 'm83_' + chnl[i] + '_' + str(index[i][j])
			filenames = [im_name, kern_name]
			choices.append(filenames) #choices-----im, kern
			k += 1
	
	##==========plot==========
	#sigma_lam, sigma_par, sigma_per = kern(wvl)
	#fig = plt.figure()
	#plt.plot(wvl, sigma_par, color='g', label='Parallel')
	#plt.plot(wvl, sigma_per, color='c', label='Perpendicular')
	#plt.plot(wvl, sigma_lam, color='k', label='Geometric mean')
	#plt.title('PSF profile')
	#plt.xlabel('Wavelength, $\lambda$ [$\mu$m]')
	#plt.ylabel('$\sigma_{\lambda}$ [MJy/sr]')
	#plt.legend(loc='upper left')
	#
	#plt.savefig(plot_save_path + 'kernel_PSF_profile.jpeg')
	#plt.show()
	
	##==========write csv==========
	with open('kern.csv', 'w') as csvfile:
	    writer = csv.writer(csvfile)
	    writer.writerow(["Images", "Kernels"])
	    writer.writerows(choices)



def conv():
	system('idl conv.pro')
