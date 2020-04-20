import numpy as np 
import math
import scipy.integrate as integrate
import histograms as hist


## ================================================================= ##
# Wasserstein distance via Sinkhorn's algorithm

def inner_prod(a, b):
	#inner product
	return np.sum(a*b)


def cum_sum(a, t):
	'''
	Function to calculate the cumulative sum of a from 0 to t.

	a - ndarray
	t - float in [0,1]
	'''

	return sum(a[:round(a.size*t)])


def inv_cum_sum(a, p):
	'''
	Function that calculates the inverse cum sum of a for a given
	probability p.

	a - ndarray
	p - float in [0,1]
	'''

	cumsum = 0; i = 0
	while cumsum < p:
		cumsum += a[i]; i += 1

	return i/a.size


def cumul_alg(a, b):
	'''
	Algorithm that calculates the Wasserstein distance between two 
	histograms a and b.
	See https://www.slideshare.net/gpeyre/an-introduction-to-optimal-transport
	(slide 48)

	a, b - ndarray of same size
	
	The algorithm calculates the Wasserstein distance by first 
	calculating the cumulative functions C_a(t) = int[-inf -> t]( da(x) )
	and the inverse cumulative functions.
	Then 
	T = C_a^-1 * C_b
	and
	W_2(a,b) = (int[0->1]( |C_a^-1 - C_b^-1| ))^(1/p)
	'''

	#Normalize a and b:
	a /= np.sum(a)
	b /= np.sum(b)

	C_a = lambda x: cum_sum(a, x)
	C_a_inv = lambda x: inv_cum_sum(a, x)
	C_b = lambda x: cum_sum(b, x)
	C_b_inv = lambda x: inv_cum_sum(b, x)

	T = lambda x: inv_cum_sum(a, cum_sum(b, x))

	W_2 = lambda x: (C_a_inv(x) - C_b_inv(x))**2


	return math.sqrt(integrate.quad(W_2, 0, 1)[0])



### SETUP
#choose two histograms a and b:
a = hist.gaussian(400, 1, 0.05, 0)
# a = hist.gaussian(400, 0.2, 0.06) + hist.gaussian(400, 0.5, 0.06)*2 + hist.gaussian(400, 0.8, 0.06)
# a = hist.slater(400, 0.5, 30, 0)
# a = hist.from_func(400, lambda x: ((x-0.5)*10)**4) + 3*hist.gaussian(400, 0.7, 0.1)
# a = hist.from_func(400, func)
# a = hist.from_func(400, lambda x: 1-x,0)
# a = hist.from_func(400, lambda x: np.cos(x*7*3.14)+1)
# a = hist.from_func(400, lambda x: x)
# a = hist.dirac_delta(400, 0.5)

b = hist.gaussian(400, 0, 0.05, 0)
# b = hist.gaussian(400, 0.8, 0.05, 0) + hist.gaussian(400, 0.2, 0.05, 0)
# b = hist.gaussian(400, 0.2, 0.06)*2 + hist.gaussian(400, 0.5, 0.06) + hist.gaussian(400, 0.8, 0.06)*2
# b = hist.from_func(400, lambda x: np.cos(x*5*3.14)+1)
# b = hist.from_func(400, lambda x: ((x-0.5)*10)**4) + 3*hist.gaussian(400, 0.7, 0.1)
# b = hist.from_func(400, func2)
# b = hist.from_func(400, lambda x: x,0)
# b = hist.from_func(400, lambda x: x**0)
# b = hist.lorentzian(400, 0.5, 0.1)


print(cumul_alg(a, b))