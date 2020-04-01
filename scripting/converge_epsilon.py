import numpy as np 
import matplotlib.pyplot as plt
import modules.histograms as hist
import modules.sinkhorn_algorithm as sink
import modules.plot as plot
import math


## ================================================================= ##
# CONVERGE_EPSILON
#script to find best epsilon to fit histograms to each other



#some histogram functions
def func(x):
	return np.where(x < 0.5, x, 1-x)

def func2(x):
	return 0.5 - func(x)

### SETUP
#choose two histograms a and b:
a = hist.gaussian(600, 0.2, 0.06)
# a = hist.gaussian(600, 0.2, 0.06) + hist.gaussian(600, 0.5, 0.06)*2 + hist.gaussian(600, 0.8, 0.06)
# a = hist.slater(600, 0.5, 30, 0)
# a = hist.from_func(400, lambda x: 1-x**2)
# a = hist.from_func(400, func)
a = hist.from_func(500, lambda x: 1-x,0)
# a = hist.from_func(600, lambda x: np.cos(x*10*3.14)+1)
# a = hist.from_func(600, lambda x: x)
# a = hist.dirac_delta(600, 0.5)

# b = hist.gaussian(400, 0.5, 0.05, 0)
b = hist.gaussian(600, 0.2, 0.06)*2 + hist.gaussian(600, 0.5, 0.06) + hist.gaussian(600, 0.8, 0.06)*2
# b = hist.from_func(600, lambda x: np.sin(x*10*3.14)+1)
# b = hist.from_func(400, func2)
# b = hist.from_func(500, lambda x: x,0)
# b = hist.from_func(600, lambda x: x**0)

