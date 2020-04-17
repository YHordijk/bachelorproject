import numpy as np 
import math
import modules.sinkhorn_algorithm as sink


## ================================================================= ##
# Barycenter calculation based on input histograms and 


def barycenter(a, weights, epsilon=0.05, cost_fn=None, error_fn=None, max_iter=100, converge_thresh=10**-5):
	'''
	Method that finds a barycenter b of hists a to solve the following:
	min_b \sum^R_{k=1} c_k W_2(a_k,b)
	given histograms a_k in hists a and weights c_k

	a - numpy matrix of histograms with histograms on rows
	weights - vector containing the weights for the minimization

	The algorithm itself is as follows:

	v0 <- np.ones(len(a[0]), len(a)) #array of R vectors.
	'''

	#default error function for two distributions x1 and x2
	if error_fn is None:
		error_fn = lambda x1, x2: np.max(abs(x1-x2))
	#default cost function for a sparse matrix specified by x and y
	if cost_fn is None:
		cost_fn = lambda x1, x2: abs(x1-x2)**2

	#normalize weights
	weights /= np.sum(weights)

	R = a.shape[0]
	N = a.shape[1]

	#normalize the hists a
	for hist in a:
		hist /= np.sum(hist)

	#initialization
	v = np.ones((R,N)).astype(float)
	u = v.copy()

	
	Y, X = np.meshgrid(np.linspace(0,1,N), np.linspace(0,1,N))
	C = cost_fn(X,Y)
	K1 = np.exp(-C/epsilon)
	K = lambda x: np.dot(np.dot(K1,x),K1)
	In = np.ones(N)

	#algorithm
	i = 0
	converged = False
	while i < max_iter and not converged:
		## =========== ##
		#iteration part one
		for k in range(R):
			u[k] = a[k]/K(v[k])

		## =========== ##
		#iteration part two
		logb = np.zeros(N).astype(float)
		for k in range(R):
			logb += weights[k] * np.log(u[k]*K(v[k]))
		b = np.exp(logb)

		## =========== ##
		#iteration part three
		for k in range(R):
			v[k] = b/K(u[k])

		## =========== ##
		#check for convergence:
		converged = 1
		P = np.empty((N,N,R))
		for k in range(R):
			P[:,:,k] = np.diag(u[k]) @ K(np.diag(v[k]))
			converged *= error_fn(a[k], P[:,:,k]@In) < converge_thresh and error_fn(b, P[:,:,k].T@In) < converge_thresh

		i += 1

	return b