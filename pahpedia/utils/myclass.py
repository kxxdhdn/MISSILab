#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

General function modes

"""

import os, sys

class blockPrint:

	def __enter__(self):
		self._original_stdout = sys.stdout
		sys.stdout = open(os.devnull, 'w')

	def __exit__(self, exc_type, exc_val, exc_tb):
		sys.stdout.close()
		sys.stdout = self._original_stdout

"""
------------------------------ MAIN (test) ------------------------------
"""
if __name__ == "__main__":

	with blockPrint():
		print("blocked print")
	print("enabled print")
