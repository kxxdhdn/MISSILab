#!/usr/bin/env python
# -*- coding: utf-8 -*-

import time, sys
t0 = time.time()

import os, sys

## TEST progressbar
##------------------
N1 = int(1E5)
N2 = 10
a = [1]*N1

from tqdm import tqdm, trange
# for i,j in tqdm(enumerate(a)):
for i in trange(N1, desc='1st loop'):
	if i%1E4==0:
		# print('coucou', i)
		tqdm.write('coucou {}'.format(i)) # tqdm does not support print()
	for k in trange(N2, desc='2nd loop', leave=False):
		if i%1E4==0:
			if k%8==0:
				tqdm.write('loulou {}'.format(k))
		# time.sleep(1)
		pass

# import progressbar
# widgets=[
# 	'1st loop: ', progressbar.Bar(), ' ', 
# 	progressbar.Percentage(), ' ', 
# 	# progressbar.Counter(), 
# 	progressbar.SimpleProgress(), ' ', 
# 	# progressbar.Timer(), ' ', 
# 	progressbar.ETA(), 
# 	# progressbar.AbsoluteETA(), 
# ]

# bar = progressbar.ProgressBar(
# 	max_value=N1*2,
# 	widgets=widgets, 
# 	redirect_stdout=True, 
# )

## progressbar does not support multi-line bars
# w2 = ['2nd loop: ', progressbar.Bar(),]
# bar2 = progressbar.ProgressBar(max_value=N2, widgets=w2, redirect_stdout=True)

# bar.start()
# for i in range(N1):
# # for i,j in bar(enumerate(a)): # don't need start(), update(), finish()
# 	if i%1E4==0:
# 		print('coucou', i)
# 	# print(i)
# 	# time.sleep(.01)
# 	bar.update(i)
# bar.finish()

## Cannot print text if flush
# import sys
# for i in range(100):
# 	time.sleep(1)
# 	# print('coucou')
# 	sys.stdout.write("\r{}%".format(i))
# 	sys.stdout.flush()

## TEST os.path
##--------------
print('Script: ', os.path.abspath(__file__))
print('Script dir: ', os.path.dirname(os.path.abspath(__file__)))
print('Current working dir: ', os.getcwd())

t_total = time.time()
print("****** total_time = {:.0f} seconds ******".format(t_total - t0))
