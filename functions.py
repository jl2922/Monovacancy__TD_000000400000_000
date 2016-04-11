# Commonly used functions
# Author: Junhao Li <streaver91@gmail.com>
# Version: 0.0.2

import sys
from scipy.optimize import fmin
from scipy.interpolate import interp1d
import numpy as np
import datetime

def partialD(fn, x0, index, dx = 1.0e-3):
	# Calculate partial derivative with respect to index at x0
	dx0 = np.zeros_like(x0)
	dx0[index] = dx
	fnP2 = fn(x0 + 2 * dx0)
	fnP = fn(x0 + dx0)
	fnM = fn(x0 - dx0)
	fnM2 = fn(x0 - 2 * dx0)
	partialD = (-1 * fnP2 + 8 * fnP - 8 * fnM + fnM2) / (12 * dx)
	return partialD

def fmin_jh(
	func,
	x0,
	args = (),
	xtol = 0.0001,
	ftol = 0.0001,
	maxfun = None,
	full_output = 0,
	disp = 1,
	retall = 0,
	callback = None
):
	def fshifted(*args):
		argsNew = (args[0] + x0,) + args[1:]
		return func(*argsNew)
	res = fmin(
		fshifted,
		np.zeros_like(x0),
		args = args,
		xtol = xtol,
		ftol = ftol,
		maxfun = maxfun,
		full_output = full_output,
		disp = disp,
		retall = retall,
		callback = callback
	)
	resNew = (res[0] + x0,) + res[1:]
	return resNew

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

def printDict(obj):
	for key in obj:
		print key + ':'
		print ' ', obj[key]

def clock(msg):
	print datetime.datetime.now().time(), msg

def _getAngle(vec1, vec2):
    # Get angle between two vectors in degrees (always between 0 - 180)
    vec1Unit = vec1 / np.linalg.norm(vec1)
    vec2Unit = vec2 / np.linalg.norm(vec2)
    angle = np.arccos(np.dot(vec1Unit, vec2Unit))
    if np.isnan(angle):
        return 0.0
    angle = angle * 180.0 / np.pi
    return angle

def rootSumSquare(*args):
	args = np.array(args)
	print args
	return np.sqrt(np.sum(args**2))
