import numpy as np 
import matplotlib.pyplot as plt


## ================================================================= ##
# plot




def plot(res, save_to=None, show_plot=True, plot_subtitle=None):
	#get results
	a, b = res.a, res.b
	P = res.P
	v, u = res.v, res.u
	error_a, error_b = res.error_a, res.error_b
	W = res.W
	bc_map = res.bc_map
	epsilon = res.epsilon
	converge_thresh = res.converge_thresh

	### PLOTTING
	fig = plt.figure(figsize=(16,9))
	title = f'Sinkhorn algorithm errors and coupling matrix for two histograms with $\epsilon={epsilon:.3f}$\nWasserstein distance $W_\epsilon={W:.2f}$; Converged after {len(error_a)} iterations with threshold$={converge_thresh}$'
	if plot_subtitle is None:
		plt.suptitle(title)
	else:
		plt.suptitle(f'{title}\n{plot_subtitle}')

	
	#plot the histograms:
	#histogram a
	plt.subplot(2,2,1)
	plt.title(f'Histogram $a$ and $P_\epsilon b$, $n={len(a)}$')
	plt.xlabel('Bin')
	plt.ylabel('Mass')
	plt.plot(range(len(a)), a, 'blue', linewidth=3, label='$a$')
	plt.plot(range(len(a)), P@b/np.sum(P@b), 'red', label='$P_\epsilon b$')
	plt.fill_between(range(len(a)), a, facecolor='lightblue')
	plt.yticks([])
	plt.xlim(0,len(a))
	plt.ylim(0,max(max(P@b/np.sum(P@b)), max(a))*1.1)
	plt.legend()

	#histogram b
	plt.subplot(2,2,2)
	plt.title(f'Histogram $b$ and $P_\epsilon^⊤ a$, $n={len(b)}$')
	plt.xlabel('Bin')
	plt.ylabel('Mass')
	plt.plot(range(len(b)), b, 'blue', linewidth=3, label='$b$')
	plt.plot(range(len(b)), P.T@a/np.sum(P.T@a), 'red', label='$P_\epsilon^⊤ a$')
	plt.fill_between(range(len(b)), b, facecolor='lightblue')
	plt.yticks([])
	plt.xlim(0,len(b))
	plt.ylim(0,max(max(P.T@a/np.sum(P.T@a)), max(b))*1.1)
	plt.legend()

	#Plot error
	plt.subplot(2,2,3)
	plt.title('Log errors during algorithm')
	plt.xlabel('Iteration')
	plt.ylabel('$log_{10}(Err)$')
	plt.plot(range(len(error_a)), np.log(error_a)/np.log(10), label='$log_{10}(Err(P\mathbb{I}  , a$))')
	plt.plot(range(len(error_b)), np.log(error_b)/np.log(10), label='$log_{10}(Err(P^⊤\mathbb{I}, b$))')
	plt.legend()

	#Plot coupling matrix
	plt.subplot(2,2,4)
	plt.title('Coupling matrix with barycentric map ($P_\epsilon$)')
	plt.xlabel('Bin in histogram b')
	plt.ylabel('Bin in histogram a')
	plt.imshow(np.log(P+1e-5), aspect='auto')

	#Plot baryocentric map:
	plt.plot(bc_map, range(len(a)), 'r', linewidth=3)
	plt.gca().invert_yaxis()

	#save plot if a file path is given
	if save_to is not None:
		plt.savefig(save_to)

	#show plot if needed
	if show_plot:
		plt.show()
	#otherwise clear and close plot to save memory
	else:
		plt.clf()
		plt.cla()
		plt.close()

	