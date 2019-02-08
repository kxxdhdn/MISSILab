import numpy as np
from scipy import optimize
from astropy.io import fits
import matplotlib.pyplot as plt


def rdata(filename):
	hdul = fits.open(filename + '.fits')
	hdr = hdul[0].header
	data = hdul[0].data
	return hdr, data

def delaberr(a, threshold=1.5):
	med_a = np.median(a)
	diff = np.abs(a - med_a)
	med_diff = np.median(diff)
	if med_diff == 0:
		s = 0
	else:
		s = diff / float(med_diff)
		#print(s)
	mask = s > threshold
	a = np.ma.masked_array(np.array(a), np.array(mask))
	#print(diff)
	#print(a)
	return a, mask

def f1(x, A, B):
	return A * x + B

def f_err(x, y, ux, uy):
	return np.sqrt(ux**2 - x * uy**2) / y


hdr, i62 = rdata('I62_val')
i77 = rdata('I77_val')[1]
i86 = rdata('I86_val')[1]
i112 = rdata('I112_val')[1]
i127 = rdata('I127_val')[1]
i170 = rdata('I170_val')[1]
i62_err = rdata('I62_unc')[1]
i77_err = rdata('I77_unc')[1]
i86_err = rdata('I86_unc')[1]
i112_err = rdata('I112_unc')[1]
i127_err = rdata('I127_unc')[1]
i170_err = rdata('I170_unc')[1]
nx = hdr['NAXIS1']
ny = hdr['NAXIS2']

##-----rapports-----
r1 = []
r2 = []
r3 = []
r4 = []
r5 = []
r1_err = []
r2_err = []
r3_err = []
r4_err = []
r5_err = []
for y in range(ny):
	for x in range(nx):
		if (i62[y,x]!=0 and i77[y,x]!=0 and i86[y,x]!=0 and i112[y,x]!=0):
			r1.append(i62[y,x] / i112[y,x])
			r2.append(i77[y,x] / i112[y,x])
			r3.append(i86[y,x] / i112[y,x])
			r4.append(i77[y,x] / i62[y,x])
			r5.append(i86[y,x] / i62[y,x])
			r1_err.append(f_err(i62[y,x], i112[y,x], i62_err[y,x], i112_err[y,x]))
			r2_err.append(f_err(i77[y,x], i112[y,x], i77_err[y,x], i112_err[y,x]))
			r3_err.append(f_err(i86[y,x], i112[y,x], i86_err[y,x], i112_err[y,x]))
			r4_err.append(f_err(i77[y,x], i62[y,x], i77_err[y,x], i62_err[y,x]))
			r5_err.append(f_err(i86[y,x], i62[y,x], i86_err[y,x], i62_err[y,x]))
#--------------------

##delete aberrant points
r1, mask = delaberr(r1)
#print(mask)
r2 = np.ma.masked_array(np.array(r2), np.array(mask))
r3 = np.ma.masked_array(np.array(r3), np.array(mask))
r4 = np.ma.masked_array(np.array(r4), np.array(mask))
r5 = np.ma.masked_array(np.array(r5), np.array(mask))
r1_err = np.ma.masked_array(np.array(r1_err), np.array(mask))
r2_err = np.ma.masked_array(np.array(r2_err), np.array(mask))
r3_err = np.ma.masked_array(np.array(r3_err), np.array(mask))
r4_err = np.ma.masked_array(np.array(r4_err), np.array(mask))
r5_err = np.ma.masked_array(np.array(r5_err), np.array(mask))

##-----plot-----
x = r3
y = r2
plt.figure(figsize=(7,5))
plt.errorbar(x, y, xerr=r3_err, yerr=r2_err, fmt='.', ms=0.8, c='k', ecolor='k')

##linear fit
popt, pcov = optimize.curve_fit(f1, x, y)
print(popt)
perr = np.sqrt(np.diag(pcov))
print(perr)
x1 = np.arange(0.1, 1, 0.0001)
#plt.plot(x1, f1(x1, *popt), c='gray')
plt.errorbar(x1, f1(x1, *popt), xerr=None, yerr=perr[1], c='white', lw=0.5, ecolor='gray')

##legend
plt.title('PAH bands ratio correlation of M83')
plt.xscale('log')
plt.yscale('log')
plt.xlim((0.08,2))
plt.ylim((0.7,10))
plt.xlabel('$I_{8.6}$/$I_{11.3}$')
plt.ylabel('$I_{7.7}$/$I_{11.3}$')

plt.savefig('correlation.jpeg')
plt.show()

