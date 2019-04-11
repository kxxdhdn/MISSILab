#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Read & write .h5 file

"""

import numpy as np
import h5py as H5

path = "./"
#path = "/Users/abricot/Github/"


def write_hdf5(filename, *dataset):

	hf = H5.File(path+filename+'.h5', 'w')
	for dset in dataset:
		hf.create_dataset(dset[0] ,data=dset[1])

	hf.flush()
	hf.close()

def read_hdf5(filename, *dset_name):

	hf = H5.File(path+filename+'.h5', 'r')
	dset_data = []
	for name in dset_name:
		data = hf.get(name)
		flag = 1
		if flag==1:
			data = np.array(data)
		dset_data.append(data)

	return dset_data

if __name__ == "__main__":

	## in situ test
	write_hdf5('test', ['dset1', [1,0,2]], ['dset2', [2,2,2]])
	data = read_hdf5('test', 'dset2', 'dset1')
	print(data)