# Commonly used functions
# Author: Junhao Li <streaver91@gmail.com>

import sys
from scipy.optimize import fmin
from scipy.interpolate import interp1d
import numpy as np
import datetime

import config as C

def fmin_jh(func, x0, args = (), xtol = 0.0001, ftol = 0.0001):
	def fshifted(x):
		return func(x + x0)
	res = fmin(
		fshifted,
		np.zeros_like(x0),
		args = args,
		xtol = xtol,
		ftol = ftol
	)
	res[0] += x0
	return res

def fmax_jh(func, x0, args = (), xtol = 0.0001, ftol = 0.0001):
	def fshifted(x):
		return -func(x + x0)
	res = fmin(
		fshifted,
		np.zeros_like(x0),
		args = args,
		xtol = xtol,
		ftol = ftol
	)
	res[0] += x0
	return res

def clock(msg):
	if C.DEBUG:
		print datetime.datetime.now().time(), msg
