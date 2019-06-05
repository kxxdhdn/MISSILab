#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from rwfits import *

source = 'm82'
chnl = 'ch2'
aorkey = '18376448'
Nexpo = 9
Ndce = 1
bcd_path = '/Users/dhu/data/spitzer/'+source+'/r'+aorkey+'/'+chnl+'/bcd/'

# (w, v)
x, y = 50, 108
flux = []
for j in range(Ndce):
	for i in range(Nexpo):
		bcd_name = 'SPITZER_S2_'+aorkey+'_'+\
		'0'*(4-len(str(i)))+str(i)+'_'+\
		'0'*(4-len(str(j)))+str(j)+'_7_bcd'
		flux.append(read_fits(bcd_path+bcd_name, None)[0][y,x])

time = np.arange(Nexpo*Ndce)
fig, ax = plt.subplots()
ax.plot(time, flux)
plt.show()
