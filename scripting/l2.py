import modules.ir as ir
import modules.histograms as hist
import modules.plot as plot
import modules.sinkhorn_algorithm as sink
import modules.colour_maps as cmap
import os
import numpy as np
import modules.comp_funcs as comp_funcs
import matplotlib.pyplot as plt 
import moviepy.editor as mvp
import ot, math
import scipy.optimize as opt




a = hist.gaussian(500, 0, 0.0001, xlim=(0,1))
a = a/np.sum(a)
distl2 = []
distws = []
for d in np.linspace(0,1,100):
	b = hist.gaussian(500, d, 0.0001, xlim=(0,1))
	b = b/np.sum(b)
	distl2.append(math.sqrt(np.sum((a-b)**2)))
	distws.append(sink.sinkhorn_test(a,b).W)
	# dist.append(np.linalg.norm(a,b))

# plt.suptitle('L2 distance')
plt.plot(np.linspace(0,1,100), distl2, label='$L_2-distance$')
plt.plot(np.linspace(0,1,100), distws, label='$W_2-distance$')
plt.xlabel('$F_a - F_b$')
plt.ylabel('D(a,b)')
plt.legend()
plt.show()