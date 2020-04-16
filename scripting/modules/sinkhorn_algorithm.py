import numpy as np 
import math


## ================================================================= ##
# Wasserstein distance via Sinkhorn's algorithm


def inner_prod(a,b):
	#inner product
	return np.sum(a*b)

def entropy(a):
	#entropy calculations
	return -np.sum(a*(np.log(a)-1))




class Results:
	def __init__(self, a, b, K, v, u, P, W, error_a, error_b, epsilon, converge_thresh):
		self.a = a 
		self.b = b
		self.K = K
		self.v = v
		self.u = u
		self.P = P
		self.W = W
		self.error_a = error_a
		self.error_b = error_b
		self.epsilon = epsilon
		self.converge_thresh = converge_thresh

		#calculate barycentric map
		self.bc_map = np.dot(K, v*np.arange(len(b))) * u / a





def sinkhorn(a, b, epsilon=0.4, cost_fn=None, error_fn=None, max_iter=100000, converge_thresh=10**-10):
	'''
	Function that implements the Sinkhorn algorithm to obtain the optimal coupling matrix P
	between two distributions a and b (in this case normalized histograms).
	We first calculate the ground-cost-matrix C followed by the Gibbs kernel K.

	The algorithm is as follows:
	v0 = np.ones(len(a))
	u <- a/np.dot(K, v)
	v <- b/np.dot(K.T, u)

	We can check our progress by calculating the error between the current approximate 
	distributions a_app and b_app and our original a and b. We obtain a_app and b_app as 
	follows:
	a_app = u * (K@v)
	b_app = v * (K.T@u)

	Numpy uses * for elementwise mult and @ for matrix mult.
	
	The error is given by np.max(abs(a-a_app)).


	a, b 			- histograms one dimensional arrays of sizes n and m
	epsilon 		- value of epsilon to use, determines strength of entropic regularisation
	cost_fn 		- function used to calculate ground cost matrix, must accept two arguments,
					  nxm and mxn arrays. Defaults to squared euclidean distance.
	error_fn		- function used to calculate error between two histograms of size n, must 
					  accept two one dimensional arrays of size n. Defaults to error.
	max_iter		- maximum number of iterations allowed for the algorithm
	converge_thresh - threshold used to determine convergence of algorithm.
					  If error < converge_thresh it is considered converged.
	'''


	#default error function for two distributions x1 and x2
	if error_fn is None:
		error_fn = lambda x1, x2: np.max(abs(x1-x2))
	#default cost function for a sparse matrix specified by x and y
	if cost_fn is None:
		cost_fn = lambda x1, x2: abs(x1-x2)**2

	## =========== ##
	#make sure a and b are normalized:
	a /= np.sum(a)
	b /= np.sum(b)


	## =========== ##
	#calculate ground-cost-matrix
	#we use the function specified in the args to calculate the cost matrix
	Y, X = np.meshgrid(np.linspace(0,1,len(a)), np.linspace(0,1,len(b)))
	C = cost_fn(X,Y)

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
		# P = np.diag(u) @ K @ np.diag(v)
		
		approx_b = v * (np.dot(K.T, u))
		error_b.append(error_fn(b, approx_b))

		## =========== ##
		#iteration part two
		v = b/(K.T @ u)
		#calculate error
		approx_a = u * (np.dot(K, v))
		error_a.append(error_fn(a, approx_a))

		#check if the algorithm has converged
		converged = error_a[-1] < converge_thresh and error_b[-1] < converge_thresh

		i += 1


	P = np.diag(u) @ K @ np.diag(v)
	W = math.sqrt(inner_prod(C,P))

	return Results(a, b, K, v, u, P, W, error_a, error_b, epsilon, converge_thresh)