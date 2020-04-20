import numpy as np
import modules.histograms as hist
import modules.sinkhorn_algorithm as sink
import modules.plot as plot
import modules.barycenter as bc


## ================================================================= ##
#MAIN
#main function used for plotting and performing sinkhorn algorithm




#histogram functions
def func(x):
	return np.where(x < 0.5, x, 1-x)

def func2(x):
	return 0.5 - func(x)

### SETUP
#choose two histograms a and b:
a = hist.gaussian(400, 0.8, 0.03, 0)
# a = hist.gaussian(400, 0.2, 0.06) + hist.gaussian(400, 0.5, 0.06)*2 + hist.gaussian(400, 0.8, 0.06)
# a = hist.slater(400, 0.5, 30, 0)
# a = hist.from_func(400, lambda x: ((x-0.5)*10)**4) + 3*hist.gaussian(400, 0.7, 0.1)
# a = hist.from_func(400, func)
# a = hist.from_func(400, lambda x: 1-x,0)
# a = hist.from_func(400, lambda x: np.cos(x*7*3.14)+1)
# a = hist.from_func(400, lambda x: x)
# a = hist.dirac_delta(400, 0.5)

b = hist.gaussian(400, 0.2, 0.06, 0)
# b = hist.gaussian(400, 0.8, 0.05, 0) + hist.gaussian(400, 0.2, 0.05, 0)
# b = hist.gaussian(400, 0.2, 0.06)*2 + hist.gaussian(400, 0.5, 0.06) + hist.gaussian(400, 0.8, 0.06)*2
# b = hist.from_func(400, lambda x: np.cos(x*5*3.14)+1)
# b = hist.from_func(400, lambda x: ((x-0.5)*10)**4) + 3*hist.gaussian(400, 0.7, 0.1)
# b = hist.from_func(400, func2)
# b = hist.from_func(400, lambda x: x,0)
# b = hist.from_func(400, lambda x: x**0)
# b = hist.lorentzian(400, 0.5, 0.1)



cost_fn = lambda x1, x2: abs(x1-x2)**10 #function used for calculating ground cost matrix
error_fn = lambda x1, x2: np.sqrt(np.sum((x1-x2)**2)) #function used for calculating errors
epsilon = 0.0001
save_to = None
plot_subtitle = None
converge_thresh = 10**-10



# barycenter = bc.barycenter(np.vstack((a,b)), (0,1))

# plot.plot_hists((a,b,barycenter), labels=('a', 'b', 'barycenter'))


### CALCULATION
#get coupling matrix and errors
res = sink.sinkhorn(a, b, epsilon, converge_thresh=converge_thresh, cost_fn=cost_fn, error_fn=error_fn)

# bc_map = res.bc_map

# plot.plot_transport(a,b,bc_map, 0.5)

plot.plot_results(res, save_to=save_to, plot_subtitle=plot_subtitle)