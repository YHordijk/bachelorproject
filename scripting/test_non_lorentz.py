import modules.ir as ir
import modules.plot as plot
import modules.sinkhorn_algorithm as sink
import os
import numpy as np
import ot
import cv2
import scipy.spatial.distance as dist
from scipy.stats import chisquare, entropy
import matplotlib.pyplot as plt



#RANDOM COLLECTION OF MOLECULES
kf_dft = [r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\DFT\AI1_DFT.t21",
		  r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\DFT\butane_DFT.t21",
		  r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\DFT\lactic acid.t21",
		  r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\DFT\l-alanine.t21",
		  r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\DFT\water.t21"]

kf_dftb = [r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\DFTB\AI1_DFTB.rkf",
		   r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\DFTB\butane_DFTB.rkf",
		   r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\DFTB\lactic acid.rkf",
		   r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\DFTB\l-alanine.rkf",
		   r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\DFTB\water.rkf"]



#HYDROCARBONS
# kf_dft = [r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\hydrocarbons\DFT\methane_DFT.t21",
# 		  r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\hydrocarbons\DFT\ethane_DFT.t21",
# 		  r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\hydrocarbons\DFT\propane_DFT.t21",
# 		  r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\hydrocarbons\DFT\butane_DFT.t21",
# 		  r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\hydrocarbons\DFT\pentane_DFT.t21",
# 		  r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\hydrocarbons\DFT\hexane_DFT.t21",
# 		  r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\hydrocarbons\DFT\heptane_DFT.t21",
# 		  r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\hydrocarbons\DFT\octane_DFT.t21"]

# kf_dftb = [r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\hydrocarbons\DFTB\methane_DFTB.rkf",
# 		   r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\hydrocarbons\DFTB\ethane_DFTB.rkf",
# 		   r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\hydrocarbons\DFTB\propane_DFTB.rkf",
# 		   r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\hydrocarbons\DFTB\butane_DFTB.rkf",
# 		   r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\hydrocarbons\DFTB\pentane_DFTB.rkf",
# 		   r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\hydrocarbons\DFTB\hexane_DFTB.rkf",
# 		   r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\hydrocarbons\DFTB\heptane_DFTB.rkf",
# 		   r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\hydrocarbons\DFTB\octane_DFTB.rkf"]


#GEOMETRY OPTIMISED AMINO-INDAN CONFORMERS
# kf_dft = [r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\conformers\AI1_DFT.t21",
# 		  r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\conformers\AI2_DFT.t21",
# 		  r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\conformers\AI3_DFT.t21",
# 		  r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\conformers\AI4_DFT.t21",
# 		  r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\conformers\AI5_DFT.t21",
# 		  r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\conformers\AI6_DFT.t21"]	

# kf_dftb = [r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\conformers\AI1_DFTB.rkf",
# 		   r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\conformers\AI2_DFTB.rkf",
# 		   r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\conformers\AI3_DFTB.rkf",
# 		   r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\conformers\AI4_DFTB.rkf",
# 		   r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\conformers\AI5_DFTB.rkf",
# 		   r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\conformers\AI6_DFTB.rkf"]


#SINGLEPOINT AMINO-INDAN CONFORMERS
# kf_dft = [r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\conf_unopt\AI1_DFT.t21",
# 		  r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\conf_unopt\AI2_DFT.t21",
# 		  r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\conf_unopt\AI3_DFT.t21",
# 		  r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\conf_unopt\AI4_DFT.t21",
# 		  r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\conf_unopt\AI5_DFT.t21",
# 		  r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\conf_unopt\AI6_DFT.t21"]

# kf_dftb = [r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\conf_unopt\AI1_DFTB.rkf",
# 		   r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\conf_unopt\AI2_DFTB.rkf",
# 		   r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\conf_unopt\AI3_DFTB.rkf",
# 		   r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\conf_unopt\AI4_DFTB.rkf",
# 		   r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\conf_unopt\AI5_DFTB.rkf",
# 		   r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\conf_unopt\AI6_DFTB.rkf"]



#get intensities and frequencies
freq_dft = [ir.get_freqs_intens(f)[0] for f in kf_dft]
ints_dft = [ir.get_freqs_intens(f)[1] for f in kf_dft]

freq_dftb = [ir.get_freqs_intens(f)[0] for f in kf_dftb]
ints_dftb = [ir.get_freqs_intens(f)[1] for f in kf_dftb]




print('Case (Unbalanced Wasserstein)')
d = np.zeros((len(freq_dft), len(freq_dftb)))
for i in range(len(freq_dft)):
	for j in range(len(freq_dftb)):
		fa = np.asarray(freq_dft[i])
		fb = np.asarray(freq_dftb[j])

		fa = fa/np.sum(fa)
		fb = fb/np.sum(fb)

		Y, X = np.meshgrid(fa, fb)
		C = (abs(X-Y)**2).T

		ia = np.asarray(ints_dft[i])
		ib = np.asarray(ints_dftb[j])

		# ia = ia/np.sum(ia)
		# ib = ib/np.sum(ib)

		P = ot.unbalanced.sinkhorn_unbalanced(ia, ib, C, 0.005, 10**-2)
		d[i,j] = np.sum(P*C)


		# plt.subplot(1,4,1)
		# plt.plot(np.asarray(freq_dft[i]), ia)
		# plt.subplot(1,4,2)
		# plt.plot(np.asarray(freq_dftb[j]), ib)
		# plt.subplot(1,4,3)
		# plt.imshow(C)
		# plt.subplot(1,4,4)
		# plt.imshow(P)
		# plt.show()

[print(*tuple(x)) for x in d]



print('Case (Balanced Wasserstein)')
d = np.zeros((len(freq_dft), len(freq_dftb)))
for i in range(len(freq_dft)):
	for j in range(len(freq_dftb)):
		fa = np.asarray(freq_dft[i])
		fb = np.asarray(freq_dftb[j])

		fa = fa/np.sum(fa)
		fb = fb/np.sum(fb)

		Y, X = np.meshgrid(fa, fb)
		C = (abs(X-Y)**2).T

		ia = np.asarray(ints_dft[i])
		ib = np.asarray(ints_dftb[j])

		ia = ia/np.sum(ia)
		ib = ib/np.sum(ib)

		P = ot.bregman.sinkhorn(ia, ib, C, 0.005)
		d[i,j] = np.sum(P*C)


		# plt.subplot(1,4,1)
		# plt.plot(np.asarray(freq_dft[i]), ia)
		# plt.subplot(1,4,2)
		# plt.plot(np.asarray(freq_dftb[j]), ib)
		# plt.subplot(1,4,3)
		# plt.imshow(C)
		# plt.subplot(1,4,4)
		# plt.imshow(P)
		# plt.show()

[print(*tuple(x)) for x in d]



#cosine similarity
def cossim(a,b):
	return np.dot(a,b)/(np.sum(np.abs(a))*np.sum(np.abs(b)))

print('Case (cosine similarity)')
d = np.zeros((len(freq_dft), len(freq_dftb)))
for i in range(len(freq_dft)):
	for j in range(len(freq_dftb)):
		fa = np.asarray(freq_dft[i])
		fb = np.asarray(freq_dftb[j])

		fa = fa/np.sum(fa)
		fb = fb/np.sum(fb)

		Y, X = np.meshgrid(fa, fb)
		C = (abs(X-Y)**2).T

		ia = np.asarray(ints_dft[i])
		ib = np.asarray(ints_dftb[j])

		ia = ia/np.sum(ia)
		ib = ib/np.sum(ib)

		d[i,j] = cossim(ia, ib)


		# plt.subplot(1,4,1)
		# plt.plot(np.asarray(freq_dft[i]), ia)
		# plt.subplot(1,4,2)
		# plt.plot(np.asarray(freq_dftb[j]), ib)
		# plt.subplot(1,4,3)
		# plt.imshow(C)
		# plt.subplot(1,4,4)
		# plt.imshow(P)
		# plt.show()

[print(*tuple(x)) for x in d]

