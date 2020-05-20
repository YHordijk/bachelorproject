
import numpy as np
import ot
import cv2
import scipy.spatial.distance as dist
from scipy.stats import chisquare, entropy



def wasserstein_distance(ir_dft, ir_dftb):
	#normalize spectra
	ir_dft = [i/np.sum(i) for i in ir_dft]
	ir_dftb = [i/np.sum(i) for i in ir_dftb]

	#wasserstein distances
	Y, X = np.meshgrid(np.linspace(0,1,ir_dft[0].size), np.linspace(0,1,ir_dft[0].size))
	C = abs(Y-X)**2

	# print('Case (Balanced Wasserstein)')
	d = np.zeros((len(ir_dft), len(ir_dftb)))
	for i, a in enumerate(ir_dft):
		for j, b in enumerate(ir_dftb):
			# d[i,j] = np.sum(ot.bregman.sinkhorn(a,b,C, 0.1)*C)
			# d[i,j] = np.sum(ot.unbalanced.sinkhorn_stabilized_unbalanced(a,b,C,0.0005, 10**-1)*C)
			d[i,j] = np.sum(ot.emd(a,b,C)*C)
			# d[i,j] = sink.sinkhorn(a,b, 0.0005).W
	return d

def wasserstein_distance_unbalanced(ir_dft, ir_dftb):
	Y, X = np.meshgrid(np.linspace(0,1,ir_dft[0].size), np.linspace(0,1,ir_dft[0].size))
	C = abs(Y-X)**2

	# print('Case (Unbalanced Wasserstein, reg_m=10**3)')
	d = np.zeros((len(ir_dft), len(ir_dftb)))
	for i, a in enumerate(ir_dft):
		for j, b in enumerate(ir_dftb):
			# d[i,j] = np.sum(ot.bregman.sinkhorn(a,b,C, 0.0005)*C)
			d[i,j] = np.sum(ot.unbalanced.sinkhorn_unbalanced(a,b,C,0.004, 10**2)*C)
			# res = sink.sinkhorn(a,b, 0.001)
			# d[i,j] = res.W
			# plot.plot_sink_results(res, title=f'DFT {i+1} | DFTB {j+1}')
	return d


def l2(ir_dft, ir_dftb):
	#L2-norm
	l2norm = lambda a, b: np.sqrt(np.sum((a*b)**2))

	# print('Case (L2-norm)')
	d = np.zeros((len(ir_dft), len(ir_dftb)))
	for i, a in enumerate(ir_dft):
		for j, b in enumerate(ir_dftb):
			d[i,j] = l2norm(a,b)
	return d



def diagonality(ir_dft, ir_dftb):
	#normalize spectra
	ir_dft = [i/np.sum(i) for i in ir_dft]
	ir_dftb = [i/np.sum(i) for i in ir_dftb]

	#diagonality of P https://math.stackexchange.com/questions/1392491/measure-of-how-much-diagonal-a-matrix-is
	Y, X = np.meshgrid(np.linspace(0,1,ir_dft[0].size), np.linspace(0,1,ir_dft[0].size))
	C = abs(Y-X)**2

	def dist(P):
		j = np.ones(P.shape[0])
		r = np.arange(P.shape[0])
		r2 = r**2

		n = j@P@j.T
		sum_x = r@P@j.T
		sum_y = j@P@r.T
		sum_x2 = r2@P@j.T
		sum_y2 = j@P@r2.T 
		sum_xy = r@P@r.T

		return (n*sum_xy - sum_x*sum_y)/(np.sqrt(n*sum_x2 - sum_x**2) * np.sqrt(n*sum_y2 - sum_y**2))

	# print('Case (Diagonality)')
	d = np.zeros((len(ir_dft), len(ir_dftb)))
	for i, a in enumerate(ir_dft):
		for j, b in enumerate(ir_dftb):
			# P = sink.sinkhorn(a,b, 0.003).P
			P = ot.emd(a,b,C)
			d[i,j] = dist(P)
	
	return d


def diagonality_unbalanced(ir_dft, ir_dftb):
	#diagonality of P https://math.stackexchange.com/questions/1392491/measure-of-how-much-diagonal-a-matrix-is
	Y, X = np.meshgrid(np.linspace(0,1,ir_dft[0].size), np.linspace(0,1,ir_dft[0].size))
	C = abs(Y-X)**2

	def dist(P):
		j = np.ones(P.shape[0])
		r = np.arange(P.shape[0])
		r2 = r**2

		n = j@P@j.T
		sum_x = r@P@j.T
		sum_y = j@P@r.T
		sum_x2 = r2@P@j.T
		sum_y2 = j@P@r2.T 
		sum_xy = r@P@r.T

		return (n*sum_xy - sum_x*sum_y)/(np.sqrt(n*sum_x2 - sum_x**2) * np.sqrt(n*sum_y2 - sum_y**2))

	# print('Case (Diagonality)')
	d = np.zeros((len(ir_dft), len(ir_dftb)))
	for i, a in enumerate(ir_dft):
		for j, b in enumerate(ir_dftb):
			# P = sink.sinkhorn(a,b, 0.003).P
			P = ot.unbalanced.sinkhorn_unbalanced(a,b,C,0.004, 10**2)
			d[i,j] = dist(P)
	
	return d




def bhattacharyya(ir_dft, ir_dftb):
	#normalize spectra
	ir_dft = [i/np.sum(i) for i in ir_dft]
	ir_dftb = [i/np.sum(i) for i in ir_dftb]
	# #Bhattacharyya distance
	# print('Case (Bhattacharyya)')
	d = np.zeros((len(ir_dft), len(ir_dftb)))
	for i, a in enumerate(ir_dft):
		for j, b in enumerate(ir_dftb):
			BC = np.sum(np.sqrt(a*b))
			d[i,j] = -np.log(BC)
	return d



def correlation(ir_dft, ir_dftb):
	#normalize spectra
	ir_dft = [i/np.sum(i) for i in ir_dft]
	ir_dftb = [i/np.sum(i) for i in ir_dftb]
	
	# print('Case (correlation)')
	d = np.zeros((len(ir_dft), len(ir_dftb)))
	for i, a in enumerate(ir_dft):
		for j, b in enumerate(ir_dftb):
			d[i,j] = dist.correlation(a,b)
	return d



def chi_square(ir_dft, ir_dftb):
	#normalize spectra
	ir_dft = [i/np.sum(i) for i in ir_dft]
	ir_dftb = [i/np.sum(i) for i in ir_dftb]

	# print('Case (chisquare)')
	d = np.zeros((len(ir_dft), len(ir_dftb)))

	for i, a in enumerate(ir_dft):
		for j, b in enumerate(ir_dftb):
			d[i,j] = chisquare(a, b)[0]
	return d



def kl_divergence(ir_dft, ir_dftb):
	#normalize spectra
	ir_dft = [i/np.sum(i) for i in ir_dft]
	ir_dftb = [i/np.sum(i) for i in ir_dftb]

	# print('Case (KL divergence)')
	d = np.zeros((len(ir_dft), len(ir_dftb)))
	for i, a in enumerate(ir_dft):
		for j, b in enumerate(ir_dftb):
			d[i,j] = 0.5*entropy(a, b) + 0.5*entropy(b, a)
	return d

