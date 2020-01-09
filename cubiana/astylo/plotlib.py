#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

matplotlib applications

"""
from astropy import units as u
import numpy as np
from scipy import optimize
import matplotlib as mpl
import matplotlib.pyplot as plt

from lib import f_lin, f_lin0
from plot import plotool

cmap = mpl.cm.viridis
norm = mpl.colors.Normalize(vmin=0, vmax=1)
COList = ['k', 'r', 'g', 'b', 'm', 'y', 'c']
LSList = ['-', '--', '-.', ':']

class spectrum(plotool):
	"""docstring for spectrum"""
	def __init__(self, figsize):
		super().__init__()
		
		self.figure(figsize=figsize)
		

