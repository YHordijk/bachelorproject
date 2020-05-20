import numpy as np 
import math
import scipy.integrate as integrate
import histograms as hist
import plot
import ir
import sinkhorn_algorithm as sink
import ot
import matplotlib.pyplot as plt


## ================================================================= ##
# Wasserstein distance via Sinkhorn's algorithm

def inner_prod(a, b):
	#inner product
	return np.sum(a*b)


def cum_sum(a, t):
	'''
	Function to calculate the cumulative sum of a from 0 to t.

	a - ndarray of size N
	t - float in [0,N]
	'''

	return sum(a[:t])


def inv_cum_sum(a, p):
	'''
	Function that calculates the inverse cum sum of a for a given
	probability p.

	a - ndarray of size N
	p - float in [0,N]
	'''

	cumsum = 0; i = 0
	while cumsum <= p and i >= a.size:
		cumsum += a[i]; i += 1

	return i


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

	W_2_f = lambda x: (C_a_inv(x) - C_b_inv(x))**2
	W_2 = math.sqrt(integrate.quad(W_2_f, 0, 1)[0])

	return W_2



def interpolate(a, b, t):
	'''
	Function that interpolates between a and b. The extent of interpolation
	is determined by t (t=0 => a, t=1 => b).
	'''

	# t = min(1, max(0, t)) #ensure t in [0,1]

	a /= np.sum(a)
	b /= np.sum(b)

	C_a = lambda x: cum_sum(a, x)
	C_a_inv = lambda x: inv_cum_sum(a, x)
	C_b = lambda x: cum_sum(b, x)
	C_b_inv = lambda x: inv_cum_sum(b, x)

	T = lambda x: inv_cum_sum(a, cum_sum(b, x))

	c = np.zeros(a.size)
	for i in range(c.size):
		l = (1-t)*i + t*T(i)
		k = l - math.floor(l)
		c[math.floor(l)] += (1-k)*a[i]
		c[math.floor(l)+1] += k*a[i]

	# c = ((1-t)*np.eye(c.size) + t*T())

	# return np.asarray(list(map(T, np.arange(0,a.size))))
	return c




### SETUP
#choose two histograms a and b:
# a = hist.gaussian(400, 0.8, 0.03, 0)
# a = hist.gaussian(400, 0.2, 0.06) + hist.gaussian(400, 0.5, 0.06)*2 + hist.gaussian(400, 0.8, 0.06)
# a = hist.slater(400, 0.1, 30, 0)
# a = hist.from_func(400, lambda x: ((x-0.5)*10)**4) + 3*hist.gaussian(400, 0.7, 0.1)
# a = hist.from_func(400, func)
# a = hist.from_func(400, lambda x: 1-x,0)
# a = hist.from_func(400, lambda x: np.cos(x*7*3.14)+1)
# a = hist.from_func(400, lambda x: x)
# a = hist.dirac_delta(400, 0.0)
a = ir.get_spectrum_from_kf(r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\DFT\lactic acid.t21", xlim=(0,2000), n=400)

# b = hist.gaussian(400, 0.2, 0.03, 0)
# b = hist.gaussian(400, 0.8, 0.05, 0) + hist.gaussian(400, 0.2, 0.05, 0)
# b = hist.gaussian(400, 0.2, 0.06)*2 + hist.gaussian(400, 0.5, 0.06) + hist.gaussian(400, 0.8, 0.06)*2
# b = hist.from_func(400, lambda x: np.cos(x*5*3.14)+1)
# b = hist.from_func(400, lambda x: ((x-0.5)*10)**4) + 3*hist.gaussian(400, 0.7, 0.1)
# b = hist.from_func(400, func2)
# b = hist.from_func(400, lambda x: x,0)
# b = hist.from_func(400, lambda x: x**0)
# b = hist.lorentzian(400, 0.5, 0.1)
# b = hist.slater(400, 0.9, 30, 0)
# b = hist.dirac_delta(400, 1)
b = ir.get_spectrum_from_kf(r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\DFTB\lactic acid.rkf", xlim=(0,2000), n=400)

plot.compare_sinkhorn(a, b, title='Comparison between Sinkhorn implementations \n IR spectrum of lactic acid')

# plot.plot_hists((a,b), ['a', 'b'], scatter=False)


# print(cumul_alg(a, b))

# c = interpolate(a, b, 0.5)
# plot.plot_hists((a, b), ('DFT', 'DFTB'))