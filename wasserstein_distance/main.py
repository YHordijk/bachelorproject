import numpy as np 
import matplotlib.pyplot as plt
import modules.histograms as hist
import modules.sinkhorn_algorithm as sink



## ================================================================= ##
# MAIN

#histogram functions
def func(x):
	return np.where(x < 0.5, x, 1-x)

def func2(x):
	return 0.5 - func(x)


### SETUP
#choose two histograms a and b:
# a = hist.gaussian(600, 0.2, 0.06)
# a = hist.gaussian(600, 0.2, 0.06) + hist.gaussian(600, 0.5, 0.06)*2 + hist.gaussian(600, 0.8, 0.06)
# a = hist.slater(600, 0.5, 5000,0)
# a = hist.from_func(400, lambda x: 1-x**2)
# a = hist.from_func(400, func)
# a = hist.from_func(500, lambda x: 1-x,0)
a = hist.from_func(600, lambda x: np.cos(x*10*3.14)+1)
# a = hist.from_func(600, lambda x: x)
# a = hist.dirac_delta(600, 0.5)

# b = hist.gaussian(400, 0.5, 0.05)
b = hist.gaussian(600, 0.2, 0.06)*2 + hist.gaussian(600, 0.5, 0.06) + hist.gaussian(600, 0.8, 0.06)*2
# b = hist.from_func(600, lambda x: np.sin(x*10*3.14)+1)
# b = hist.from_func(400, func2)
# b = hist.from_func(500, lambda x: x,0)
# b = hist.from_func(600, lambda x: x**0)


converge_thresh = 10**-8
epsilon = 100


### CALCULATION
#get coupling matrix and errors
W, P, error_a, error_b, bc_map = sink.sinkhorn(a, b, epsilon, converge_thresh=converge_thresh)


### PLOTTING
plt.suptitle(f'Sinkhorn algorithm errors and coupling matrix for two histograms with $\epsilon={epsilon}$\nWasserstein distance $W_\epsilon={round(W,3)}$ \nConverged after {len(error_a)} iterations with threshold$={converge_thresh}$')

#plot the histograms:
#histogram a
plt.subplot(2,2,1)
plt.title(f'Histogram a, $n={len(a)}$')
plt.xlabel('Bin')
plt.ylabel('Mass')
plt.bar(range(len(a)), a, width=1)
# plt.xticks([])
plt.yticks([])
plt.xlim(0,len(a))

#histogram b
plt.subplot(2,2,2)
plt.title(f'Histogram b, $n={len(b)}$')
plt.xlabel('Bin')
plt.ylabel('Mass')
plt.bar(range(len(b)), b, width=1)
# plt.xticks([])
plt.yticks([])
plt.xlim(0,len(b))

#Plot error
plt.subplot(2,2,3)
plt.title('Log errors during algorithm')
plt.xlabel('Iteration')
plt.ylabel('log(Err)')
plt.plot(range(len(error_a)), np.log(error_a)/np.log(10), label='$log_{10}(Err(P\mathbb{I}  , a$))')
plt.plot(range(len(error_b)), np.log(error_b)/np.log(10), label='$log_{10}(Err(P^⊤\mathbb{I}, b$))')
plt.legend()

#Plot coupling matrix
plt.subplot(2,2,4)
plt.title('Coupling matrix with barycentric map ($P_\epsilon$)')
plt.xlabel('Bin in histogram b')
plt.ylabel('Bin in histogram a')
plt.imshow(np.log(P+1e-5), aspect='auto')
#Calculate baryocentric map:
plt.plot(bc_map, range(len(a)), 'r', linewidth=3)
plt.gca().invert_yaxis()
# plt.xticks([])
# plt.yticks([])

plt.show()