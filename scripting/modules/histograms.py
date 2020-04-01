import numpy as np



## ================================================================= ##
# HISTOGRAM FUNCTIONS
#collection of function to generate histograms
#the function from_func allows you to generate a histogram based
#on any supplied function




def gaussian(n, x0, sigma, min_mass=0.02):
	'''
	Function that returns a histogram based on a gaussian function
	exp(-(bins-x0)**2/(2*sigma**2))

	n - number of bins
	x0 - center of gauss curve
	sigma - standard deviation
	min_mass - mass to add to all bins before normalization
	'''

	f = lambda bins: np.exp(-(bins-x0)**2/(2*sigma**2))
	return from_func(n, f, min_mass)



def lorentzian(n, x0, w, min_mass=0):
	'''
	Function that returns a histogram based on a lorentzian function
	(x0-x)/(w/2), used for spectrum generation.

	n - number of bins
	x0 - center of lorentzian
	w - full width at half maximum of peak
	'''

	f = lambda bins: 1/(1+((x0-bins)/(w/2))**2)
	return from_func(n, f, min_mass)



def slater(n, x0, chi, min_mass=0.02):
	'''
	Function that returns a histogram based on a slater function
	exp(-abs(bins-x0)*chi)

	n - number of bins
	x0 - center of gauss curve
	chi - exponential factor
	min_mass - mass to add to all bins before normalization
	'''

	f = lambda bins: np.exp(-abs(bins-x0)*chi)
	return from_func(n, f, min_mass)



def dirac_delta(n, x0):
	'''
	Function that returns a histogram where only one bin is filled

	n - number of bins
	x0 - position of filled bin
	'''

	bins = np.linspace(0,1,n)
	index = np.argsort(np.abs(bins-x0))[0]
	h = np.zeros(n)
	h[index] = 1
	return h
	


def from_func(n, func, min_mass=0.02):
	'''
	Function that returns a histogram based on any function

	n - number of bins
	func - any function used to calculate the histogram
	min_mass - mass to add to all bins before normalization
	'''

	#get the x-values of the bins first
	bins = np.linspace(0,1,n) 
	#apply function to get histogram
	h = func(bins)
	#add min_mass
	h += min_mass
	#normalize and return histogram mass
	return h/np.sum(h)