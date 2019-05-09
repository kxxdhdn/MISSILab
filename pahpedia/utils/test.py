import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
"""
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(8.5, 5))

for ax in axes.flat:
    ax.set_axis_off()
    im = ax.imshow(np.random.random((16, 16)), cmap='viridis',
                   vmin=0, vmax=1)


fig.subplots_adjust(bottom=0.2, top=0.7, left=0.1, right=0.85,
                    wspace=0.02, hspace=0.02)
# notice that here we use ax param of figure.colorbar method instead of

# the cax param as the above example

cbar = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.95)

cbar.set_ticks(np.arange(0, 1.1, 0.5))
cbar.set_ticklabels(['low', 'medium', 'high'])

plt.show()
"""
hdr = fits.open('../test_examples/n66_LL1_cube.fits')[0].header
for kw in hdr.keys():
	print(kw)
print("------------")
if hdr['NAXIS']==3:
	for kw in hdr.keys():
		print(kw)
		# print(kw, ('3' in kw))
