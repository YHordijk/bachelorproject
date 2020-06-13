import numpy as np
import modules.histograms as hist
import modules.sinkhorn_algorithm as sink
import modules.brute_force as bf
import modules.plot as plot
import modules.barycenter as bc
import modules.ir as ir
from scipy.stats import entropy
from time import perf_counter
import matplotlib.pyplot as plt


## ================================================================= ##
#MAIN
#main function used for plotting and performing sinkhorn algorithm




#histogram functions
def func(x):
	return np.where(x < 0.5, x, 1-x)

def func2(x):
	return 0.5 - func(x)

def func3(x):
	return np.random.rand(x.size)

### SETUP
#choose two histograms a and b:
a = 0.5*hist.gaussian(20, 0., 0.15, 0) + 0.5*hist.gaussian(20, 1, 0.15, 0)
# a = hist.from_func(400, func3)
# a = hist.gaussian(400, 0.2, 0.06) + hist.gaussian(400, 0.5, 0.06)*2 + hist.gaussian(400, 0.8, 0.06)
# a = hist.slater(400, 0.5, 30, 0)
# a = hist.from_func(400, lambda x: ((x-0.5)*10)**4) + 3*hist.gaussian(400, 0.7, 0.1)
# a = hist.from_func(400, func)
# a = hist.from_func(400, lambda x: 1-x,0)
# a = hist.from_func(400, lambda x: np.cos(x*7*3.14)+1)
# a = hist.from_func(400, lambda x: x)
# a = hist.dirac_delta(400, 0.5)
# a = ir.get_spectrum_from_kf(r"D:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\functionals\Aminoindan_CONF_2_LDA_DFT.t21", xlim=(0,4000), n=600)
# a = ir.get_spectrum_from_kf(r"D:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\functionals\Aminoindan_CONF_4_LDA_DFT.t21", xlim=(0,4000), n=10)

b = hist.gaussian(20, 0.5, 0.15, 0)
# b = hist.from_func(400, func3)
# b = hist.gaussian(400, 0.8, 0.05, 0) + hist.gaussian(400, 0.2, 0.05, 0)
# b = hist.gaussian(400, 0.2, 0.06)*2 + hist.gaussian(400, 0.5, 0.06) + hist.gaussian(400, 0.8, 0.06)*2
# b = hist.from_func(400, lambda x: np.cos(x*5*3.14)+1)
# b = hist.from_func(400, lambda x: ((x-0.5)*10)**4) + 3*hist.gaussian(400, 0.7, 0.1)
# b = hist.from_func(400, func2)
# b = hist.from_func(400, lambda x: x,0)
# b = hist.from_func(400, lambda x: x**0)
# b = hist.lorentzian(400, 0.5, 0.1)
# b = ir.get_spectrum_from_kf(r"D:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\functionals\Aminoindan_CONF_2_OLYP_DFT.t21", xlim=(0,4000), n=600)
# b = ir.get_spectrum_from_kf(r"D:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\functionals\Aminoindan_CONF_4_DFTB3_freq_DFTB.rkf", xlim=(0,4000), n=10)



cost_fn = lambda x1, x2: abs(x1-x2)**2 #function used for calculating ground cost matrix
error_fn = lambda x1, x2: np.sqrt(np.sum((x1-x2)**2)) #function used for calculating errors
epsilon = 0.001
save_to = None
plot_subtitle = None
converge_thresh = 10**-6







# vals = []
# thresh_range = np.linspace(0,8, 9)
# for i in thresh_range:
# 	print(i)
# 	vals.append(bf.wasserstein_bf(a,b, error_thresh=10**-i)[1])
bf.wasserstein_bf(a,b, error_thresh=10**-6)[1]

# plt.plot(thresh_range, vals)
# plt.show()
res = sink.sinkhorn(a, b, epsilon, converge_thresh=converge_thresh, max_iter=100000)
plt.imshow(res.P)
plt.show()
print(res.W)









# plot.plot_hist(ir.get_spectrum_from_kf(r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\second_set\Aminoindan_CONF_1_DFT.t21", xlim=(0,4000), n=600), invert_xaxis=True)

# barycenter = bc.barycenter(np.vstack((a,b)), (0,1))

# plot.plot_hists((a,b,barycenter), labels=('a', 'b', 'barycenter'))


### CALCULATION
# #get coupling matrix and errors
# for n in [300, 600, 900, 1200, 1500, 1800, 2100, 2400]:
# 	a = ir.get_spectrum_from_kf(r"D:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\functionals\Aminoindan_CONF_2_LDA_DFT.t21", xlim=(0,4000), n=n)
# 	b = ir.get_spectrum_from_kf(r"D:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\functionals\Aminoindan_CONF_2_OLYP_DFT.t21", xlim=(0,4000), n=n)

# 	start = perf_counter()
# 	for i in range(5):
# 		res = sink.sinkhorn(a, b, epsilon, converge_thresh=converge_thresh, max_iter=100000)
# 	print(f'Sinkhorn: n={n}; average time={(perf_counter()-start)/5}')


# for n in [300, 600, 900, 1200, 1500, 1800, 2100, 2400]:
# 	a = ir.get_spectrum_from_kf(r"D:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\functionals\Aminoindan_CONF_2_LDA_DFT.t21", xlim=(0,4000), n=n)
# 	b = ir.get_spectrum_from_kf(r"D:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\functionals\Aminoindan_CONF_2_OLYP_DFT.t21", xlim=(0,4000), n=n)

# 	start = perf_counter()
# 	for i in range(5):
# 		res = sink.sinkhorn_test(a, b, epsilon, converge_thresh=converge_thresh, max_iter=100000)
# 	print(f'Exact: n={n}; average time={(perf_counter()-start)/5}')

# bc_map = res.bc_map

# plot.plot_transport(a,b,bc_map, 0.5)
# res = sink.sinkhorn(a, b, epsilon, converge_thresh=converge_thresh, max_iter=100000)
# print(entropy(a,b)/2 + entropy(b,a)/2)
# plot.plot_sink_results(res, save_to=save_to, plot_subtitle=plot_subtitle)

# import math
# import matplotlib.pyplot as plt
# import matplotlib.transforms as trns
# import matplotlib.gridspec as gridspec

# l2 = []
# W2 = []



# d = 0.8
# a = hist.gaussian(400, 0, 0.03, 0, xlim=(-.3,1.3))
# b = hist.gaussian(400, d, 0.03, 0, xlim=(-.3,1.3))
# a /= np.sum(a)
# b /= np.sum(b)
# plot.plot_hists((a,b), labels=('a', 'b'))

# for d in np.linspace(0,0.8,100):
# 	a = hist.gaussian(400, 0, 0.03, 0, xlim=(-0.3,1.3))
# 	b = hist.gaussian(400, d, 0.03, 0, xlim=(-0.3,1.3))

# 	a /= np.sum(a)
# 	b /= np.sum(b)

# 	print(np.sum(np.abs(a-b)))
# 	l2.append((np.sum((np.abs(a-b))**2)))
# 	W2.append(sink.sinkhorn(a, b, 10**-2).W)

# print(math.sqrt(np.sum(a)**2+np.sum(b)**2))

# plt.plot(np.linspace(0,0.8,100), l2, label='L2-distance')
# plt.plot(np.linspace(0,0.8,100), W2, label='W2-distance')
# plt.xlabel('$m_b - m_a$')
# plt.ylabel('Distance metric')
# plt.legend()
# plt.show()



# gridspec.GridSpec(4,4)

# plt.subplot2grid((4,4), (0,1), colspan=3, rowspan=1)
# plt.plot(np.linspace(0,1,400), b)
# plt.xlabel = ''
# plt.ylabel = ''
# plt.xticks([])
# plt.yticks([])

# plt.subplot2grid((4,4), (1,0), colspan=1, rowspan=3)
# base = plt.gca().transData
# rot = trns.Affine2D().rotate_deg(90)
# plt.plot(np.linspace(0,1,400), a, transform= rot + base)
# plt.xlabel = ''
# plt.ylabel = ''
# plt.xticks([])
# plt.yticks([])

# plt.subplot2grid((4,4), (1,1), colspan=3, rowspan=3)

# plt.suptitle('Optimal coupling matrix P')
# plt.xticks([])
# plt.yticks([])
# plt.xlabel = ''
# plt.ylabel = ''



# P = sink.sinkhorn(a, b, 10**-3).P
# plt.imshow(np.log(P+1e-5), aspect='auto')
# plt.gca().invert_yaxis()

# plt.show()