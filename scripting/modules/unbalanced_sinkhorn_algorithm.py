import numpy as np 
import math
from math import log


## ================================================================= ##
# Wasserstein distance via Sinkhorn's algorithm


def inner_prod(a,b):
    #inner product
    return np.sum(a*b)

def entropy(a):
    #entropy calculations
    return -np.sum(a*(np.log(a)-1))




class Results:
    def __init__(self, a, b, K, v, u, P, W, error_a, error_b, epsilon, converge_thresh, C):
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
        self.C = C

        #calculate barycentric map
        self.bc_map = np.dot(K, v*np.arange(len(b))) * u / a



def sinkhorn(a, b, epsilon=0.4, tau=0.4, cost_fn=None, error_fn=None, max_iter=100000, converge_thresh=10**-10, print_complexity=True):
    '''
    Function that implements the Sinkhorn algorithm to obtain the optimal coupling matrix P
    between two distributions a and b (in this case un-normalized histograms).
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


    a, b            - histograms one dimensional arrays of sizes n and m
    epsilon         - value of epsilon to use, determines strength of entropic regularisation
    cost_fn         - function used to calculate ground cost matrix, must accept two arguments,
                      nxm and mxn arrays. Defaults to squared euclidean distance.
    error_fn        - function used to calculate error between two histograms of size n, must 
                      accept two one dimensional arrays of size n. Defaults to error.
    max_iter        - maximum number of iterations allowed for the algorithm
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
    #calculate ground-cost-matrix
    #we use the function specified in the args to calculate the cost matrix
    Y, X = np.meshgrid(np.linspace(0,1,len(a)), np.linspace(0,1,len(b)))
    C = cost_fn(X,Y)

    


    ## =========== ##
    #setup
    sum_a = np.sum(a)
    sum_b = np.sum(b)

    a_log = np.log(a)
    b_log = np.log(b)

    assert a.size == b.size

    n = a.size

    S = (sum_a + sum_b)/2 + 0.5 + 1/(4*log(n))
    T = (sum_a + sum_b)/2 * (log((sum_a + sum_b)/2) + 2*log(n) - 1) + log(n) + 2.5

    U = max(S+T, 2*epsilon, 4*epsilon*log(n)/tau, 4*epsilon*(sum_a + sum_b)*log(n)/tau)
    eta = epsilon/U

    A = C/eta
    A = A - A.min(axis=0)
    A = (A.T - A.min(axis=1)).T

    #calculate the gibbs kernel:
    K = np.exp(-A)

    scale_factor = (eta*tau) / (eta+tau)

    #initialize f and g as identity vectors
    f = np.zeros_like(a)
    g = np.zeros_like(b)

    P = np.diag(np.exp(f/eta)) @ K @ np.diag(np.exp(g/eta))


    ## =========== ##
    #algorithm part

    #follow deviations during iterations
    error_a = []
    error_b = []

    converged = False
    i = 0
    while i < max_iter and not converged:
        f_prev = f
        g_prev = g

        approx_a = P.sum(axis=1)
        f = (f/eta + a_log - np.log(approx_a)) * scale_factor
        P = np.diag(np.exp(f/eta)) @ K @ np.diag(np.exp(g/eta))

        approx_b = P.sum(axis=0)
        g = (g/eta + b_log - np.log(approx_b)) * scale_factor
        P = np.diag(np.exp(f/eta)) @ K @ np.diag(np.exp(g/eta))


        error_b.append(error_fn(b, approx_b))
        error_a.append(error_fn(a, approx_a))

        #check if the algorithm has converged
        converged = error_a[-1] < converge_thresh and error_b[-1] < converge_thresh
        # converged = error_a[-1] < converge_thresh

        i += 1


    W = math.sqrt(inner_prod(C,P))

    return Results(a, b, K, f, g, P, W, error_a, error_b, epsilon, converge_thresh, C)