#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Read & write .csv file

"""

import csv
import numpy as np

def read_csv(filename):

	with open(filename+'.csv', 'r') as csvfile:
		reader = csv.reader(csvfile)
		data = []
		for line in reader:
			data.append(line)

	return np.array(data)

def write_csv(filename, tags, lists):

	with open(filename+'.csv', 'w') as csvfile:
		writer = csv.writer(csvfile)
		writer.writerow(tags)
		writer.writerows(lists)

"""
------------------------------ MAIN (test) ------------------------------
"""
if __name__ == "__main__":

	## in situ test
	path = '../test_examples/'
	filename = 'write_test'

	a = np.arange(2)
	a = [a]*10
	a = np.array(a)
	print("a = ", a.shape)
	write_csv(path+filename, ['a0', 'a1'], a)
	data = read_csv(path+filename)
	print(data.shape)
