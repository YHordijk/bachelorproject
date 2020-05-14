import modules.ir as ir
import modules.plot as plot
import modules.sinkhorn_algorithm as sink
import os
import numpy as np
import ot
import cv2
import scipy.spatial.distance as dist
from scipy.stats import chisquare, entropy


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



# #HYDROCARBONS
# kf_dft = [r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\[21-05](13-05-2020)\hydrocarbons_methane_DFT.t21",
# 		  r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\[21-05](13-05-2020)\hydrocarbons_ethane_DFT.t21",
# 		  r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\[21-05](13-05-2020)\hydrocarbons_propane_DFT.t21",
# 		  r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\[21-05](13-05-2020)\hydrocarbons_butane_DFT.t21",
# 		  r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\[21-05](13-05-2020)\hydrocarbons_pentane_DFT.t21",
# 		  r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\[21-05](13-05-2020)\hydrocarbons_hexane_DFT.t21",
# 		  r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\[21-05](13-05-2020)\hydrocarbons_heptane_DFT.t21",
# 		  r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\[21-05](13-05-2020)\hydrocarbons_octane_DFT.t21"]

# kf_dftb = [r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\[21-05](13-05-2020)\hydrocarbons_methane_DFTB.rkf",
# 		   r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\[21-05](13-05-2020)\hydrocarbons_ethane_DFTB.rkf",
# 		   r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\[21-05](13-05-2020)\hydrocarbons_propane_DFTB.rkf",
# 		   r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\[21-05](13-05-2020)\hydrocarbons_butane_DFTB.rkf",
# 		   r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\[21-05](13-05-2020)\hydrocarbons_pentane_DFTB.rkf",
# 		   r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\[21-05](13-05-2020)\hydrocarbons_hexane_DFTB.rkf",
# 		   r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\[21-05](13-05-2020)\hydrocarbons_heptane_DFTB.rkf",
# 		   r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\[21-05](13-05-2020)\hydrocarbons_octane_DFTB.rkf"]


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





#get the spectra
ir_dft = [ir.get_spectrum_from_kf(f,width=1) for f in kf_dft]
ir_dftb = [ir.get_spectrum_from_kf(f,width=1) for f in kf_dftb]

# #normalize spectra
# ir_dft = [i/np.sum(i) for i in ir_dft]
# ir_dftb = [i/np.sum(i) for i in ir_dftb]


#wasserstein distances
Y, X = np.meshgrid(np.linspace(0,1,600), np.linspace(0,1,600))
C = abs(Y-X)**2


# print('Case (Balanced Wasserstein)')
# d = np.zeros((len(ir_dft), len(ir_dftb)))
# for i, a in enumerate(ir_dft/np.sum(ir_dft)):
# 	for j, b in enumerate(ir_dftb/np.sum(ir_dftb)):
# 		# d[i,j] = np.sum(ot.bregman.sinkhorn(a,b,C, 0.1)*C)
# 		# d[i,j] = np.sum(ot.unbalanced.sinkhorn_stabilized_unbalanced(a,b,C,0.0005, 10**-1)*C)
# 		d[i,j] = np.sum(ot.emd(a,b,C)*C)
# 		# d[i,j] = sink.sinkhorn(a,b, 0.0005).W
# [print(*tuple(x)) for x in d]


print('Case (Unbalanced Wasserstein, reg_m=10**3)')
d = np.zeros((len(ir_dft), len(ir_dftb)))
for i, a in enumerate(ir_dft):
	for j, b in enumerate(ir_dftb):
		# d[i,j] = np.sum(ot.bregman.sinkhorn(a,b,C, 0.0005)*C)
		d[i,j] = np.sum(ot.unbalanced.sinkhorn_unbalanced(a,b,C,0.004, 10**2)*C)
		# d[i,j] = sink.sinkhorn(a,b, 0.0005).W
[print(*tuple(x)) for x in d]


# #L2-norm
# l2 = lambda a, b: np.sqrt(np.sum((a*b)**2))

# print('Case (L2-norm)')
# d = np.zeros((len(ir_dft), len(ir_dftb)))
# for i, a in enumerate(ir_dft):
# 	for j, b in enumerate(ir_dftb):
# 		d[i,j] = l2(a,b)
# [print(*tuple(x)) for x in d]




# #diagonality of P https://math.stackexchange.com/questions/1392491/measure-of-how-much-diagonal-a-matrix-is
# def diagonality(P):
# 	j = np.ones(P.shape[0])
# 	r = np.arange(P.shape[0])
# 	r2 = r**2

# 	n = j@P@j.T
# 	sum_x = r@P@j.T
# 	sum_y = j@P@r.T
# 	sum_x2 = r2@P@j.T
# 	sum_y2 = j@P@r2.T 
# 	sum_xy = r@P@r.T

# 	return (n*sum_xy - sum_x*sum_y)/(np.sqrt(n*sum_x2 - sum_x**2) * np.sqrt(n*sum_y2 - sum_y**2))


# print('Case (Diagonality)')
# d = np.zeros((len(ir_dft), len(ir_dftb)))
# for i, a in enumerate(ir_dft):
# 	for j, b in enumerate(ir_dftb):
# 		# P = sink.sinkhorn(a,b, 0.003).P
# 		P = ot.emd(a,b,C)
# 		d[i,j] = diagonality(P)
# [print(*tuple(x)) for x in d]




# # #Bhattacharyya distance
# print('Case (Bhattacharyya)')
# d = np.zeros((len(ir_dft), len(ir_dftb)))
# for i, a in enumerate(ir_dft):
# 	for j, b in enumerate(ir_dftb):
# 		BC = np.sum(np.sqrt(a*b))
# 		d[i,j] = -np.log(BC)
# [print(*tuple(x)) for x in d]




# print('Case (correlation)')
# d = np.zeros((len(ir_dft), len(ir_dftb)))
# for i, a in enumerate(ir_dft):
# 	for j, b in enumerate(ir_dftb):
# 		d[i,j] = dist.correlation(a,b)
# [print(*tuple(x)) for x in d]




# print('Case (chisquare)')
# d = np.zeros((len(ir_dft), len(ir_dftb)))
# for i, a in enumerate(ir_dft/np.sum(ir_dft)):
# 	for j, b in enumerate(ir_dftb/np.sum(ir_dftb)):
# 		# d[i,j] = np.sum(ot.bregman.sinkhorn(a,b,C, 0.1)*C)
# 		d[i,j] = chisquare(a, b)[0]
# 		# d[i,j] = sink.sinkhorn(a,b, 0.0005).W
# [print(*tuple(x)) for x in d]



# print('Case (KL divergence)')
# d = np.zeros((len(ir_dft), len(ir_dftb)))
# for i, a in enumerate(ir_dft/np.sum(ir_dft)):
# 	for j, b in enumerate(ir_dftb/np.sum(ir_dftb)):
# 		d[i,j] = 0.5*entropy(a, b) + 0.5**entropy(b, a)
# [print(*tuple(x)) for x in d]