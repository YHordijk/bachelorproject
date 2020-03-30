import numpy as np 
import matplotlib.pyplot as plt

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


## ================================================================= ##
# Wasserstein distance via Sinkhorn's algorithm

def sinkhorn(a, b, epsilon, max_iter=10000, converge_thresh=10**-5):
	'''
	Function that implements the Sinkhorn algorithm to obtain the optimal coupling matrix P
	between two distributions a and b (in this case normalized histograms).
	We first calculate the ground-cost-matrix C followed by the Gibbs kernel K.

	The algorithm is as follows:
	v0 = np.ones(len(a))
	u <- a/np.dot(K, v)
	v <- b/np.dot(K.T, u)

	We can check our progress by calculating the RMSD between the current approximate 
	distributions a_app and b_app and our original a and b. We obtain a_app and b_app as 
	follows:
	a_app = u * (K@v)
	b_app = v * (K.T@u)

	Numpy uses * for elementwise mult and @ for matrix mult.
	
	The RMSD is then simply given by np.sqrt(np.sum(a-a_app)**2).
	'''

	#RMSD function for two distributions x1 and x2
	RMSD = lambda x1, x2: np.sqrt(np.sum((x1-x2)**2))


	## =========== ##
	#calculate ground-cost-matrix
	#we use the square of the distance between the bins as the cost
	Y, X = np.meshgrid(np.linspace(0,1,len(a)), np.linspace(0,1,len(b)))
	C = (X-Y)**2

	#calculate the gibbs kernel:
	K = np.exp(-C/epsilon).T

	## =========== ##
	#algorithm part

	#initialize v as identity vector. We get u from v later
	Im = np.ones(len(a))
	In = np.ones(len(b))
	v = In

	#follow deviations during iterations
	error_a = []
	error_b = []

	converged = False
	i = 0
	while i < max_iter and not converged:
		## =========== ##
		#iteration part one
		u = a/(K @ v)
		#Calculate error
		approx_b = v * (np.dot(K.T, u))
		error_b.append(RMSD(b, approx_b))

		## =========== ##
		#iteration part two
		v = b/(K.T @ u)
		#calculate error
		approx_a = u * (np.dot(K, v))
		error_a.append(RMSD(a, approx_a))

		#check if the algorithm has converged
		converged = error_a[-1] < converge_thresh and error_b[-1] < converge_thresh

		i += 1


	#calculate P
	P = np.diag(u) @ K @ np.diag(v)

	#calculate baryocentric map

	bc_map = np.dot(K, v*np.linspace(0,1,len(a))) * u / a

	return P, error_a, error_b, bc_map

## ================================================================= ##
# MAIN

def func(x):
	return np.where(x < 0.5, x, 1-x)

def func2(x):
	return 0.5 - func(x)


#get two histograms:
# a = get_gaussian_histogram(300, 0.2, 0.06)
# b = get_gaussian_histogram(300, 0.6, 0.2)
# a = get_histogram(400, lambda x: 1-x**2)
a = get_histogram(400, func, 0)
# b = get_histogram(600, lambda x: np.sin(x*4*3.14)+1)
b = get_histogram(400, func2, 0)


#plot the histograms:
plt.subplot(2,2,1)
plt.title('Histogram a')
plt.bar(range(len(a)), a, width=1)
plt.subplot(2,2,2)
plt.title('Histogram b')
plt.bar(range(len(b)), b, width=1)

#get coupling matrix and errors
P, error_a, error_b, bc_map = sinkhorn(a, b, 0.01)

#Plot error
plt.subplot(2,2,3)
plt.title('Errors during iteration')
plt.plot(range(len(error_a)), error_a)
plt.plot(range(len(error_b)), error_b)

#Plot coupling matrix
plt.subplot(2,2,4)
plt.title('Coupling matrix')
plt.imshow(np.log(P+1e-5))

#Calculate baryocentric map:
plt.plot(range(len(a)), bc_map, 'r', linewidth=3)

plt.show()

