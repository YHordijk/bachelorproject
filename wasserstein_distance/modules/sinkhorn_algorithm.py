import numpy as np 
import math


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
	
	The error (RMSD) is then simply given by np.sqrt(np.sum(a-a_app)**2).
	'''

	#error function for two distributions x1 and x2
	error = lambda x1, x2: np.sqrt(np.sum((x1-x2)**2))


	## =========== ##
	#calculate ground-cost-matrix
	#we use the square of the distance between the bins as the cost
	Y, X = np.meshgrid(np.linspace(0,1,len(a)), np.linspace(0,1,len(b)))
	C = abs(X-Y)**2

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
		error_b.append(error(b, approx_b))

		## =========== ##
		#iteration part two
		v = b/(K.T @ u)
		#calculate error
		approx_a = u * (np.dot(K, v))
		error_a.append(error(a, approx_a))

		#check if the algorithm has converged
		converged = error_a[-1] < converge_thresh and error_b[-1] < converge_thresh

		i += 1


	#calculate P
	P = np.diag(u) @ K @ np.diag(v)

	#calculate barycentric map
	bc_map = np.dot(K, v*np.arange(len(b))) * u / a

	#compute wasserstein distance
	W = 0
	# for i in range(P.shape[0]):
	# 	for j in range(P.shape[1]):
	# 		W += P[i,j]*(math.log(P[i,j]/K[i,j])-1)

	return W, P, error_a, error_b, bc_map