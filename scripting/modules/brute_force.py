import scipy.optimize as opt 
import numpy as np
import matplotlib.pyplot as plt

error_fn = lambda x1, x2: np.max(abs(x1-x2))
cost_fn = lambda x1, x2: abs(x1-x2)**2


def wasserstein_bf(a, b, C=None, error_thresh=10**-10):
	'''
	We attempt to minimize the wasserstein distance given by
	W = np.sum(P, C) for optimal transport matrix P and cost
	matrix C.


	We have two inequality constraints:
	We must have that
	(1) error_fun(P.T@Ia, a) <= error_thresh
	and
	(2) error_fun(P@Ib, b) <= error_thresh


	for Ia = np.ones((a.size,1)) # one vector with same size as a
	and Ib = np.ones((b.size,1))

	Our bound is that all elements in 0 <= P <= inf

	We then simply initialize the matrix P as a uniform matrix 
	and let scipy.minimize optimize the matrix P
	'''


	#build cost matrix
	if C is None:
		X, Y = np.meshgrid(np.linspace(0,1,len(a)), np.linspace(0,1,len(b)))
		C = cost_fn(X,Y)


	#objective func
	def objective(P):
		temp = np.reshape(P, (a.size, b.size))
		return np.sum(temp.T*C)



	#Constraints
	#We need that PI = a and P.T I = b
	Ia = np.ones((a.size,1))
	Ib = np.ones((b.size,1))


	def constraint1(P):
		'''
		inequality constraint function is of form f(x) >= 0
		'''
		temp = np.reshape(P, (a.size, b.size))
		return error_thresh - error_fn(temp.T@Ia, a)
	con1 = {'type': 'ineq', 'fun': constraint1}

	def constraint2(P):
		temp = np.reshape(P, (a.size, b.size))
		return error_thresh - error_fn(temp@Ib, b)
	con2 = {'type': 'ineq', 'fun': constraint2}

	cons = [con1, con2]

	bounds = opt.Bounds(0, np.inf)

	P0 = np.ones((a.size, b.size)).flatten()/(np.sum(a) + np.sum(b)) #first guess

	opts = {'maxiter':1000, 'disp':True}

	minimised = opt.minimize(objective, P0, constraints=cons, options=opts, bounds=bounds) 

	return minimised.x.reshape((a.size, b.size)), minimised.fun



def unbalanced_wasserstein_bf(a, b, reg, reg_m):
	NotImplemented