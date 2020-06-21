import modules.ir as ir
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



def C_L2L2(Fa, Ia, Fb, Ib):
	return (Fa-Fb)**2 + (Ia-Ib)**2

def C_L2L1(Fa, Ia, Fb, Ib):
	return (Fa-Fb)**2 + abs(Ia-Ib)

def C_L1L1(Fa, Ia, Fb, Ib):
	return abs(Fa-Fb) + abs(Ia-Ib)

def C_prod(Fa, Ia, Fb, Ib):
	return (Fa-Fb)**2*(Ia-Ib)**2


def lp(Fa, Ia, Fb, Ib, C):
	def objective(P):
		Pt = P.reshape((Fa.size, Fb.size, 2))
		return np.sum(C*P)

	bound = opt.Bounds(0, np.inf)
	# constraint1 = lambda P: 



a = r"D:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\DFT\l-alanine.t21"
b = r"D:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\DFTB\l-alanine.rkf"

fa, ia = ir.get_freqs_intens(a)
fb, ib = ir.get_freqs_intens(b)

fmax = max(fa.max(), fb.max())
imax = max(ia.max(), ib.max())
fa = fa/fmax
fb = fb/fmax
ia = ia/imax
ib = ib/imax

Fa, Fb = np.meshgrid(fa, fb)
Ia, Ib = np.meshgrid(ia, ib)

C = C_L1L1(Fa, Ia, Fb, Ib)


plt.imshow(-np.log(C), aspect='auto')
plt.gca().invert_yaxis()
plt.show()



a = np.vstack((fa, ia))
b = np.vstack((fb, ib))

P = ot.emd(ia,ib,C)
plt.imshow(P)
plt.show()
