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



a = r"C:\Users\Yuman Hordijk\Desktop\Scripts\bachelorproject\scripting\RUNS\KFFiles\l-alanine.t21"
b = r"C:\Users\Yuman Hordijk\Desktop\Scripts\bachelorproject\scripting\RUNS\KFFiles\dftb_l_alanine.rkf"

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

C = C_prod(Fa, Ia, Fb, Ib)


# plt.imshow(-np.log(C), aspect='auto')
# plt.gca().invert_yaxis()
# plt.show()



a = np.vstack((fa, ia))
b = np.vstack((fb, ib))

# P = ot.emd(ia/ia.sum(),ib/ib.sum(),C)
res = sink.sinkhorn(ia, ib, cost_mat=C)
Pi = res.P 

# P = ot.emd(ia, ib, C)
Wi = np.sum(Pi*C)


# plt.imshow(Pi)
# plt.gca().invert_yaxis()
# plt.show()


# P = ot.emd(ia/ia.sum(),ib/ib.sum(),C)
res = sink.sinkhorn(fa, fb, cost_mat=C)
Pf = res.P 

# P = ot.emd(ia, ib, C)
Wf = np.sum(Pf*C)

# plt.imshow(Pf)
# plt.gca().invert_yaxis()
# plt.show()


print(math.sqrt(Wi), math.sqrt(Wf), math.sqrt(Wi+Wf))


plt.subplot(1,2,1)
plt.imshow(Pf+Pi)
plt.gca().invert_yaxis()

plt.subplot(1,2,2)
plt.imshow(Pf*Pi)
plt.gca().invert_yaxis()
plt.show()