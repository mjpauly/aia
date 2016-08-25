"""module for finding the brightness values for AIA images using the date
"""

import pickle
import numpy as np
import imp
from matplotlib.dates import date2num, num2date
import os
pb0r = imp.load_source('pb0r',
					   os.path.join(os.path.dirname(__file__), 'pb0r.py'))


def get_intensity(date, wavelength):
	"""takes in a datetime instance and gets the intensity of the Sun
	based on its apparent brightness on that day
	"""
	date = date2num(date)
	reg = openreg(wavelength)
	idx = find_nearest(reg[0],date)
	return reg[1][idx]


def find_nearest(array, value):
	idx = (np.abs(array-value)).argmin()
	return idx


def openreg(wave):
	with open(os.path.join(os.path.dirname(__file__),
						   wave + 'regression.pkl'), 'rb') as f:
		return pickle.load(f)
