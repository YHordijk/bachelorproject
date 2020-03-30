import numpy as np

## ================================================================= ##
# HISTOGRAM FUNCTIONS

def get_gaussian_histogram(n, x0, sigma, min_mass=0.02):
	'''
	Function that returns a histogram based on a gaussian function
	exp(-(bins-x0)**2/(2*sigma**2))

	n - number of bins
	x0 - center of gauss curve
	sigma - standard deviation
	min_mass - mass to add to all bins before normalization
	'''

	#get the x-values of the bins first
	bins = np.linspace(0,1,n) 
	#calculate the guassian histogram
	h = np.exp(-(bins-x0)**2/(2*sigma**2))
	#add min_mass
	h += min_mass
	#normalize and return histogram mass
	return h/np.sum(h)


def get_histogram(n, func, min_mass=0.02):
	'''
	Function that returns a histogram based on any function

	n - number of bins
	func - any function used to calculate the histogram
	min_mass - mass to add to all bins before normalization
	'''

	bins = np.linspace(0,1,n) 

	h = func(bins)
	h += min_mass

	return h/np.sum(h)