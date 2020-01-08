#!/usr/bin/env python
# -*- coding: utf-8 -*-

class base:
	"""docstring for base"""
	def __init__(self, a1=111, a2=222):
		print('base')
		self.a1 = a1
		self.a2 = a2
		print('base')

	def f0(self, x0):
		self.y0 = x0 + self.a1

class l1(base):
	"""docstring for l1"""
	def __init__(self, b1=8, b2=88, b3=888):
		print('\nlevel 1')
		super().__init__(b1, b2)
		self.b1 = b1
		self.b2 = b2
		self.b3 = b3
		print(self.a1)
		print(self.a2)
		print('level 1\n')

	def f1(self, x1):
		self.y1 = x1 + self.b1

		return self.y1

class l2(base):
	"""docstring for l2"""
	def __init__(self, c1=9, c2=99):
		print('\nlevel 2')
		super().__init__(c1, c2)
		self.c1 = c1
		self.c2 = c2
		print(self.a1)
		print(self.a2)
		print('level 2\n')

class l3(l1, base):
	"""docstring for l3"""
	def __init__(self, d1, d2, d3):
		print('\nlevel3')
		super().__init__(d1, d2)
		# base.__init__(self, 6, 66)
		# self.f0(100)
		# print(self.y0)
		self.f1(100)
		print(self.y1)
		print('------')
		print(self.a1) 
		print(self.a2)
		print(self.b1)
		print(self.b2)
		print(self.b3)
		# print(self.c1)
		# print(self.c2)
		print('------')
		print('level3\n')

	def f3(self, x3):
		print(self.f1(7))
		print(self.f1(x3))
		# self.__init__(6, 66, 666)

	def f33(self):
		super().__init__(4,44)
		self.f3(10)
		
		
if __name__ == "__main__":
	
	l3(5, 55, 555).f33()
