#!/usr/bin/env python
# -*- coding: utf-8 -*-

import time
t0 = time.time()

import os
print('Script: ', os.path.abspath(__file__))
print('Script dir: ', os.path.dirname(os.path.abspath(__file__)))
print('Current working dir: ', os.getcwd())

t_total = time.time()
print("****** total_time = {:.0f} seconds ******".format(t_total - t0))
